"""GPU MLIP screening -> selected periodic DFT single points -> relative-site comparison.

The four notebook models are optional adapters, not a claim of better predictions.
Only a fixed model is used per run. Missing models/errors never trigger substitution.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import re
import time

import numpy as np
from ase import Atoms
from ase.build import add_adsorbate, fcc111
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import BFGS

from molecular import save_dft_config
from periodic import PeriodicConfig

ROOT = Path(__file__).resolve().parent
METALS = ('Ag', 'Au', 'Cu', 'Pd', 'Pt')
LATTICE_CONSTANTS = dict(Ag=4.09, Au=4.08, Cu=3.61, Pd=3.89, Pt=3.92)
VALENCE = dict(H=1, Ag=11, Au=11, Cu=11, Pd=18, Pt=18)
MODEL_CHOICES = ('mace', 'sevennet', 'chgnet', 'mattersim')


def dump_json(data, path):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def structure_hash(atoms):
    data = dict(numbers=atoms.numbers.tolist(), positions=np.round(atoms.positions, 10).tolist(),
                cell=np.round(atoms.cell.array, 10).tolist(), pbc=atoms.pbc.tolist())
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def comparable_signature(atoms):
    """Only equal composition, cell and PBC have a meaningful relative site energy."""
    return (tuple(sorted(Counter(atoms.get_chemical_symbols()).items())),
            tuple(np.round(atoms.cell.array, 8).ravel()), tuple(atoms.pbc),
            structure_hash(atoms[:-1]) if atoms[-1].symbol == 'H' else None)


def write_table(rows, path, fields):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def group_minima(rows, directory, tolerance=.15):
    """Cluster H positions by minimum-image distance; fixed substrates must match."""
    from ase.geometry import find_mic
    groups = []
    for row in sorted((r for r in rows if r['status'] == 'converged'), key=lambda r: (r['energy_ev'], r['id'])):
        path = (Path(directory) / row['final_file']).resolve()
        if not path.is_relative_to(Path(directory).resolve()):
            raise ValueError('후보 경로가 탐색 폴더 외부입니다.')
        atoms = read(str(path), format='extxyz')
        if structure_hash(atoms) != row['final_sha256']:
            raise ValueError('완화 후 구조가 변경되었습니다.')
        signature = comparable_signature(atoms)
        representative = row['id']
        for group_id, group_atoms, group_signature in groups:
            if signature == group_signature:
                _, distance = find_mic(atoms.positions[-1] - group_atoms.positions[-1], atoms.cell, atoms.pbc)
                if float(distance) < tolerance:
                    representative = group_id
                    break
        if representative == row['id']:
            groups.append((row['id'], atoms, signature))
        row['equivalent_minimum_of'] = representative
    return len(groups)


def generate_candidates(directory, material='HEA', size=(3, 3, 4), seed=0, vacuum=9.):
    directory = Path(directory)
    if directory.exists() and any(directory.iterdir()):
        raise ValueError('후보 생성에는 비어 있는 새 출력 폴더를 사용하세요.')
    if material not in ('HEA', *METALS) or len(size) != 3 or min(size) < 1 or np.prod(size) > 199:
        raise ValueError('HEA / Ag / Au / Cu / Pd / Pt 및 1~199개 표면 원자를 지정하세요.')
    if size[2] < 2 or vacuum < 6:
        raise ValueError('표면은 최소 2층, 양쪽 진공은 각각 6 Å 이상으로 만드세요.')
    directory.mkdir(parents=True, exist_ok=True)
    a0 = np.mean(list(LATTICE_CONSTANTS.values())) if material == 'HEA' else LATTICE_CONSTANTS[material]
    slab = fcc111('Ag' if material == 'HEA' else material, size=size, a=a0, vacuum=vacuum)
    if material == 'HEA':
        if len(slab) < 5:
            raise ValueError('5원소 HEA에는 최소 5개 표면 원자가 필요합니다.')
        # Same deterministic counts/order as the supplied notebooks, without executing them.
        counts = [len(slab) // 5 + (i < len(slab) % 5) for i in range(5)]
        symbols = [e for e, count in zip(METALS, counts) for _ in range(count)]
        np.random.default_rng(seed).shuffle(symbols)
        slab.set_chemical_symbols(symbols)
    write(str(directory / 'clean.extxyz'), slab, format='extxyz')
    rows = []
    for site in ('ontop', 'bridge', 'fcc', 'hcp'):
        for i in range(size[0]):
            for j in range(size[1]):
                atoms = slab.copy()
                add_adsorbate(atoms, 'H', 1.8, position=site, offset=(i, j))
                atoms.set_constraint(FixAtoms(indices=list(range(len(slab)))))
                name = f'{site}_{i}_{j}'
                write(str(directory / f'{name}.extxyz'), atoms, format='extxyz')
                stored = read(str(directory / f'{name}.extxyz'), format='extxyz')
                rows.append(dict(id=name, initial_site=site, file=f'{name}.extxyz',
                                 structure_sha256=structure_hash(stored)))
    manifest = dict(schema_version=1, material=material, size=list(size), seed=seed,
                    lattice_constant_angstrom=float(a0), vacuum_each_side_angstrom=vacuum,
                    substrate_composition=dict(Counter(slab.get_chemical_symbols())),
                    relaxation='adsorbate_only', candidates=rows)
    dump_json(manifest, directory / 'candidates.json')
    return directory / 'candidates.json'


def make_calculator(model, device='cuda', checkpoint=None, acceleration='none'):
    """Explicit notebook model/version. Checkpoint override is only for MACE."""
    import torch
    if device == 'cuda' and not torch.cuda.is_available():
        raise RuntimeError('CUDA를 사용할 수 없습니다. CPU 실행을 원하면 --device cpu를 명시하세요.')
    if checkpoint and model != 'mace':
        raise ValueError('--checkpoint는 현재 MACE에서만 지원합니다.')
    if acceleration not in ('none', 'cueq') or (acceleration == 'cueq' and (model != 'mace' or device != 'cuda')):
        raise ValueError('cueq 가속은 MACE / CUDA 조합에서만 사용할 수 있습니다.')
    torch.set_num_threads(4)
    metadata = dict(adapter=model, device=device, torch_version=torch.__version__, acceleration=acceleration,
                    hardware=torch.cuda.get_device_name(0) if device == 'cuda' else 'CPU')
    started = time.perf_counter()
    if model == 'mace':
        os.environ.setdefault('XDG_CACHE_HOME', str(ROOT / '.model_cache'))
        from mace.calculators import mace_mp
        from mace.calculators.foundations_models import download_mace_mp_checkpoint
        resolved_checkpoint = Path(checkpoint) if checkpoint else Path(download_mace_mp_checkpoint('medium'))
        metadata['checkpoint_sha256'] = hashlib.sha256(resolved_checkpoint.read_bytes()).hexdigest()
        calc = mace_mp(model=str(resolved_checkpoint), default_dtype='float64', device=device, enable_cueq=acceleration == 'cueq')
        metadata.update(package_version=version('mace-torch'), model='MACE-MP medium', dtype='float64')
        if acceleration == 'cueq':
            metadata['cuequivariance_version'] = version('cuequivariance')
            metadata['cuequivariance_ops_version'] = version('cuequivariance-ops-torch-cu12')
        from mace.calculators.foundations_models import mace_mp_urls
        metadata['checkpoint_source'] = str(checkpoint) if checkpoint else mace_mp_urls['medium']
        # MACE retains model objects; fingerprint actual tensor values, not just a model alias.
        digest = hashlib.sha256()
        for key, value in sorted(calc.models[0].state_dict().items()):
            digest.update(key.encode())
            digest.update(value.detach().cpu().contiguous().numpy().tobytes())
        metadata['weights_sha256'] = digest.hexdigest()
    elif model == 'sevennet':
        from sevenn.calculator import SevenNetCalculator
        calc = SevenNetCalculator(model='7net-0', device=device)
        metadata.update(package_version=version('sevenn'), model='7net-0')
    elif model == 'chgnet':
        from chgnet.model.model import CHGNet
        from chgnet.model.dynamics import CHGNetCalculator
        calc = CHGNetCalculator(model=CHGNet.load(model_name='0.3.0'), use_device=device)
        metadata.update(package_version=version('chgnet'), model='CHGNet 0.3.0')
    elif model == 'mattersim':
        from mattersim.forcefield import MatterSimCalculator
        calc = MatterSimCalculator(load_path='MatterSim-v1.0.0-5M.pth', device=device)
        metadata.update(package_version=version('mattersim'), model='MatterSim-v1.0.0-5M.pth')
    else:
        raise ValueError(f'알 수 없는 모델: {model}')
    metadata['load_seconds'] = time.perf_counter() - started
    return calc, metadata


def select_candidates(rows, top_k, explore=1, seed=0):
    good = sorted((r for r in rows if r['status'] == 'converged' and r.get('equivalent_minimum_of', r['id']) == r['id']),
                  key=lambda r: (r['energy_ev'], r['id']))
    if not good:
        return []
    selected = good[:top_k]
    rest = good[top_k:]
    if rest and explore:
        indices = np.random.default_rng(seed).choice(len(rest), min(explore, len(rest)), replace=False)
        selected += [rest[int(i)] for i in indices]
    return [r['id'] for r in selected]


def screen_candidates(manifest_path, model='mace', device='cuda', top_k=3, explore=1,
                      fmax=.03, steps=150, seed=0, checkpoint=None, calculator=None, acceleration='none'):
    if top_k < 1 or explore < 0 or not 0 < fmax <= 1 or steps < 0:
        raise ValueError('top-k ≥ 1, explore ≥ 0, 0 < fmax ≤ 1, steps ≥ 0이 필요합니다.')
    manifest_path = Path(manifest_path).resolve()
    manifest = read_json(manifest_path)
    ids = [r['id'] for r in manifest['candidates']]
    if len(ids) != len(set(ids)) or not all(re.fullmatch(r'[A-Za-z0-9_-]{1,80}', name) for name in ids):
        raise ValueError('후보 ID는 중복 없는 영문·숫자·밑줄·하이픈이어야 합니다.')
    directory = manifest_path.parent
    output = directory / ('screen_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    output.mkdir()
    calc, metadata = calculator if calculator is not None else make_calculator(model, device, checkpoint, acceleration)
    rows, signature = [], None
    started = time.perf_counter()
    for candidate in manifest['candidates']:
        tick = time.perf_counter()
        row = dict(candidate, status='error', energy_ev=None, max_free_force_ev_a=None,
                   max_total_force_ev_a=None, steps=0, error=None)
        try:
            path = (directory / candidate['file']).resolve()
            if not path.is_relative_to(directory):
                raise ValueError('후보 파일은 manifest 폴더 내부에 있어야 합니다.')
            atoms = read(str(path), format='extxyz')
            if structure_hash(atoms) != candidate['structure_sha256']:
                raise ValueError('생성 후 입력 구조가 변경되었습니다.')
            this_signature = comparable_signature(atoms)
            if signature is not None and signature != this_signature:
                raise ValueError('조성·셀·주기 경계가 다른 후보를 같은 순위로 비교할 수 없습니다.')
            signature = this_signature
            if atoms[-1].symbol != 'H':
                raise ValueError('마지막 원자가 H인 흡착 후보가 필요합니다.')
            initial_positions = atoms.positions.copy()
            atoms.set_constraint(FixAtoms(indices=list(range(len(atoms) - 1))))
            atoms.calc = calc
            optimizer = BFGS(atoms, logfile=str(output / (candidate['id'] + '.log')), maxstep=.15)
            converged = bool(optimizer.run(fmax=fmax, steps=steps))
            energy = float(atoms.get_potential_energy())
            free_force = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
            total_force = float(np.linalg.norm(atoms.get_forces(apply_constraint=False), axis=1).max())
            if not np.isfinite([energy, free_force, total_force]).all():
                raise RuntimeError('유한하지 않은 MLIP 결과입니다.')
            if not np.allclose(atoms.positions[:-1], initial_positions[:-1], atol=1e-10):
                raise RuntimeError('고정된 표면 원자가 이동했습니다.')
            # Require the independently checked force criterion as well as ASE's return value.
            converged = converged and free_force < fmax
            atoms.calc = None
            final = output / (candidate['id'] + '.extxyz')
            write(str(final), atoms, format='extxyz')
            row.update(status='converged' if converged else 'not_converged', energy_ev=energy,
                       max_free_force_ev_a=free_force, max_total_force_ev_a=total_force,
                       steps=optimizer.nsteps, final_file=final.name, final_sha256=structure_hash(read(str(final), format='extxyz')))
        except Exception as exc:
            row['error'] = str(exc)
        row['runtime_seconds'] = time.perf_counter() - tick
        rows.append(row)
        print(f"{candidate['id']}: {row['status']} E={row['energy_ev']} · {row['runtime_seconds']:.2f}s", flush=True)
        # Every completed candidate survives interruption of a later candidate.
        dump_json(dict(model=metadata, rows=rows, status='running'), output / 'screen.json')
    minimum_count = group_minima(rows, output)
    selected = select_candidates(rows, top_k, explore, seed)
    good = [r['energy_ev'] for r in rows if r['status'] == 'converged']
    for r in rows:
        r['relative_energy_ev'] = r['energy_ev'] - min(good) if r['status'] == 'converged' else None
    result = dict(schema_version=1, status='complete' if good else 'failed', model=metadata,
                  manifest=str(manifest_path), rows=rows, selected_ids=selected,
                  selection=dict(top_k=top_k, random_exploration=explore, seed=seed, minimum_separation_angstrom=.15),
                  distinct_minima=minimum_count,
                  fmax_ev_a=fmax, max_steps=steps, relaxation='H_only_fixed_substrate',
                  elapsed_seconds=time.perf_counter() - started,
                  caveats=['MLIP 수렴은 DFT 정확도 검증이 아닙니다.',
                           'initial_site는 시작 위치이며 완화 후 다른 사이트로 이동할 수 있습니다.',
                           '탐색하지 않은 표면 배치·coverage·스핀은 순위에 포함되지 않습니다.'])
    dump_json(result, output / 'screen.json')
    write_table(rows, output / 'screen.csv', ['id', 'initial_site', 'status', 'energy_ev', 'relative_energy_ev',
                'max_free_force_ev_a', 'max_total_force_ev_a', 'steps', 'runtime_seconds', 'equivalent_minimum_of', 'error'])
    return output / 'screen.json'


def reselect_screen(screen_path, top_k=3, explore=1, seed=0):
    """Reuse completed energies/forces without repeating inference or overwriting a run."""
    screen_path = Path(screen_path)
    data = read_json(screen_path)
    if data.get('status') != 'complete' or top_k < 1 or explore < 0:
        raise ValueError('완료된 탐색과 올바른 후보 개수가 필요합니다.')
    data['distinct_minima'] = group_minima(data['rows'], screen_path.parent)
    data['selected_ids'] = select_candidates(data['rows'], top_k, explore, seed)
    data['selection'] = dict(top_k=top_k, random_exploration=explore, seed=seed, minimum_separation_angstrom=.15)
    data['reselected_from'] = screen_path.name
    output = screen_path.with_name('screen_selected_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '.json')
    dump_json(data, output)
    return output


def relative_comparison(rows):
    good = [r for r in rows if r['status'] == 'converged']
    if len(good) < 2:
        return dict(status='insufficient_converged_pairs', count=len(good))
    best = min(good, key=lambda r: r['dft_energy_ev'])
    errors = [(r['ml_energy_ev'] - best['ml_energy_ev']) - (r['dft_energy_ev'] - best['dft_energy_ev']) for r in good]
    pairs = [(a, b) for i, a in enumerate(good) for b in good[i + 1:]
             if abs(a['dft_energy_ev'] - b['dft_energy_ev']) > .01 and abs(a['ml_energy_ev'] - b['ml_energy_ev']) > .01]
    return dict(status='compared_selected_subset', count=len(good),
                dft_best_id=best['id'], ml_best_id=min(good, key=lambda r: r['ml_energy_ev'])['id'],
                relative_energy_mae_ev=float(np.mean(np.abs(errors))),
                pair_order_agreement=(sum((a['dft_energy_ev'] - b['dft_energy_ev']) * (a['ml_energy_ev'] - b['ml_energy_ev']) > 0 for a, b in pairs) / len(pairs)) if pairs else None,
                pair_tie_tolerance_ev=.01,
                scope='Selected structures only; DFT basis/kmesh convergence not established.')


def verify_selected(screen_path, backend='auto', settings=None, limit=None, progress=print):
    from dft_workflow import run_calculation, CalculationError
    screen_path = Path(screen_path).resolve()
    screen = read_json(screen_path)
    if limit is not None and limit < 1:
        raise ValueError('검증 개수는 1 이상이어야 합니다.')
    if screen.get('status') != 'complete':
        raise ValueError('완료된 탐색 결과가 필요합니다.')
    selected = screen['selected_ids'] if limit is None else screen['selected_ids'][:limit]
    if not selected:
        raise ValueError('검증할 수렴 후보가 없습니다.')
    output = screen_path.parent / ('dft_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    output.mkdir()
    settings = dict(settings or {})
    unknown = set(settings) - (set(PeriodicConfig.__dataclass_fields__) - {'atoms', 'lattice', 'system'})
    if unknown:
        raise ValueError('DFT 수치 설정에 허용되지 않는 필드: ' + ', '.join(sorted(unknown)))
    rows, signature = [], None
    for candidate_id in selected:
        if not re.fullmatch(r'[A-Za-z0-9_-]{1,80}', candidate_id):
            raise ValueError('잘못된 후보 ID입니다.')
        source = next(r for r in screen['rows'] if r['id'] == candidate_id)
        row = dict(id=candidate_id, ml_energy_ev=source['energy_ev'], status='error')
        try:
            final = (screen_path.parent / source['final_file']).resolve()
            if not final.is_relative_to(screen_path.parent):
                raise ValueError('후보 경로가 탐색 폴더 외부입니다.')
            atoms = read(str(final), format='extxyz')
            if source['status'] != 'converged' or structure_hash(atoms) != source['final_sha256']:
                raise ValueError('미수렴 또는 변경된 후보는 검증할 수 없습니다.')
            current = comparable_signature(atoms)
            if signature is not None and signature != current:
                raise ValueError('서로 다른 조성·셀·경계 조건을 비교할 수 없습니다.')
            signature = current
            config_settings = dict(settings)
            if 'spin' not in config_settings:
                config_settings['spin'] = sum(VALENCE[e] for e in atoms.get_chemical_symbols()) % 2
            config_settings.setdefault('kmesh', (2, 2, 1))
            config = PeriodicConfig.from_ase(atoms, **config_settings)
            result, directory = run_calculation(config, output / candidate_id, backend, progress)
            row.update(status=result['status'], dft_energy_ev=result['extrapolated_energy_ev'],
                       result_directory=str(directory), dft_seconds=result['elapsed_seconds'],
                       config=config.to_dict())
        except (ValueError, RuntimeError, OSError, KeyError) as exc:
            row['error'] = str(exc)
            if isinstance(exc, CalculationError):
                row['result_directory'] = str(exc.directory)
        rows.append(row)
        dump_json(dict(rows=rows, comparison=relative_comparison(rows)), output / 'comparison.json')
    comparison = dict(rows=rows, comparison=relative_comparison(rows), energy='(E + F) / 2, eV per cell',
                      source_screen=str(screen_path),
                      note='MLIP 완화 구조에서 DFT 단일점 검증. 절대 흡착에너지가 아니며 H2 기준은 포함하지 않음.')
    dump_json(comparison, output / 'comparison.json')
    write_table(rows, output / 'comparison.csv', ['id', 'status', 'ml_energy_ev', 'dft_energy_ev', 'dft_seconds', 'result_directory', 'error'])
    return output / 'comparison.json'


def main():
    parser = argparse.ArgumentParser(description='표면 H 후보 생성 → MLIP 탐색 → 선택 후보 DFT 검증')
    sub = parser.add_subparsers(dest='command', required=True)
    gen = sub.add_parser('generate')
    gen.add_argument('--output', type=Path, required=True)
    gen.add_argument('--material', choices=['HEA', *METALS], default='HEA')
    gen.add_argument('--size', nargs=3, type=int, default=(3, 3, 4))
    gen.add_argument('--seed', type=int, default=0)
    screen = sub.add_parser('screen')
    screen.add_argument('manifest', type=Path)
    screen.add_argument('--model', choices=MODEL_CHOICES, default='mace')
    screen.add_argument('--device', choices=['cuda', 'cpu'], default='cuda')
    screen.add_argument('--checkpoint', type=Path)
    screen.add_argument('--acceleration', choices=['none', 'cueq'], default='none')
    screen.add_argument('--top-k', type=int, default=3)
    screen.add_argument('--explore', type=int, default=1)
    screen.add_argument('--fmax', type=float, default=.03)
    screen.add_argument('--steps', type=int, default=150)
    screen.add_argument('--seed', type=int, default=0)
    verify = sub.add_parser('verify')
    verify.add_argument('screen', type=Path)
    verify.add_argument('--backend', choices=['auto', 'local', 'wsl'], default='auto')
    verify.add_argument('--settings', type=Path, help='수치 설정 JSON: kmesh, grid_level, spin 등 (구조 제외)')
    verify.add_argument('--limit', type=int)
    select = sub.add_parser('select', help='완료된 탐색의 중복 위치를 제거하고 후보 재선정')
    select.add_argument('screen', type=Path)
    select.add_argument('--top-k', type=int, default=3)
    select.add_argument('--explore', type=int, default=1)
    select.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    try:
        if args.command == 'generate':
            output = generate_candidates(args.output, args.material, tuple(args.size), args.seed)
        elif args.command == 'screen':
            output = screen_candidates(args.manifest, args.model, args.device, args.top_k, args.explore,
                                       args.fmax, args.steps, args.seed, args.checkpoint, acceleration=args.acceleration)
        elif args.command == 'select':
            output = reselect_screen(args.screen, args.top_k, args.explore, args.seed)
        else:
            output = verify_selected(args.screen, args.backend, read_json(args.settings) if args.settings else None, args.limit)
        print(f'결과: {output}', flush=True)
        data = read_json(output)
        return 2 if data.get('status') == 'failed' or any(r['status'] != 'converged' for r in data.get('rows', [])) else 0
    except (ValueError, RuntimeError, ImportError, OSError) as exc:
        print(f'실행 실패: {exc}', flush=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
