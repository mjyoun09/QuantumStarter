"""Geometry optimization: geomeTRIC (TRIC internal coordinates) proposes steps,
this project's engines supply the energy and analytic gradient at each geometry.

L.-P. Wang, C. Song, J. Chem. Phys. 144, 214108 (2016). Every step is a full
recorded calculation in steps/step_NNN; the final geometry is recomputed with the
requested output settings in the run directory. A step whose SCF does not
converge stops the optimization: its forces would be meaningless.
"""
from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import numpy as np

from molecular import OPT_CONVERGENCE, Atom, parse_scan, save_dft_config
from scf_control import check_cancelled
from .pyscf_engine import BOHR_ANGSTROM, HARTREE_TO_EV, calculate_molecule


class StepNotConverged(RuntimeError):
    pass


def _xyz(atoms, comment):
    return f'{len(atoms)}\n{comment}\n' + ''.join(
        f'{a.element:<2} {a.position[0]:.10f} {a.position[1]:.10f} {a.position[2]:.10f}\n' for a in atoms)


def _recording_engine(base, molecule, step_config, directory, progress, label):
    """geomeTRIC engine whose energy and gradient come from this project's engines.

    Every call is a full recorded run in steps/step_NNN, started from the previous
    call's density. Returns (engine, trajectory rows, state: last density/atoms, gradient per step)."""
    initial = step_config.atoms
    steps_dir = directory / 'steps'
    trajectory, state = [], {'density': None, 'atoms': initial, 'gradients': {}}

    def forward(message):
        check_cancelled(directory)  # The GUI writes cancel.request into the run directory.
        if message.startswith(('계산 오류', 'SCF 미수렴')):
            progress(message)

    class Engine(base):
        def calc_new(self, coords, dirname):
            check_cancelled(directory)
            step = len(trajectory) + 1
            positions = np.asarray(coords, dtype=float).reshape(-1, 3) * BOHR_ANGSTROM
            atoms = tuple(Atom(a.element, tuple(p)) for a, p in zip(initial, positions))
            step_dir = steps_dir / f'step_{step:03d}'
            summary = calculate_molecule(replace(step_config, atoms=atoms), step_dir, forward, state['density'])
            if not summary['converged'] or summary['gradient'] is None:
                raise StepNotConverged(f'{step}단계 SCF가 수렴하지 않아 {label} 계산을 중단합니다. '
                                       f'{step_dir}의 진단 결과를 확인하세요.')
            with np.load(step_dir / 'electronic_structure.npz', allow_pickle=False) as data:
                state['density'] = np.asarray(data['density_matrix'])
            state['atoms'] = atoms
            gradient = np.asarray(summary['gradient']['gradient_hartree_per_bohr'])
            state['gradients'][step] = gradient.ravel()
            energy = summary['total_energy_hartree']
            row = {'step': step, 'energy_hartree': energy,
                   'energy_change_hartree': energy - trajectory[-1]['energy_hartree'] if trajectory else None,
                   'rms_gradient_hartree_per_bohr': float(np.sqrt(np.mean(gradient**2))),
                   'max_gradient_hartree_per_bohr': float(np.abs(gradient).max()),
                   'max_atomic_force_ev_per_angstrom': summary['gradient']['max_atomic_force_ev_per_angstrom'],
                   'scf_cycles': summary['scf_cycles'], 'directory': f'steps/step_{step:03d}',
                   # geomeTRIC's own label: Hessian displacements (transition-state search) vs. steps.
                   'geometric_dirname': str(dirname), 'kind': 'hessian' if 'hess' in str(dirname).lower() else 'step',
                   'coordinates_angstrom': positions.tolist()}
            trajectory.append(row)
            change = '—' if row['energy_change_hartree'] is None else f'{row["energy_change_hartree"]:+.3e}'
            progress(f'{label} {step:3d}단계 · E = {energy:.10f} Ha · ΔE {change} · '
                     f'최대 힘 {row["max_atomic_force_ev_per_angstrom"]:.5f} eV/Å · SCF {row["scf_cycles"]}회')
            progress('OPT_JSON ' + json.dumps({k: v for k, v in row.items() if k != 'coordinates_angstrom'}))
            return {'energy': energy, 'gradient': gradient.ravel()}

    return Engine(molecule), trajectory, state


def _start(config, directory, purpose):
    """Shared start of geomeTRIC jobs: run folder, request record, Å geometry, per-step settings."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / 'restart_input.npz').exists():
        raise ValueError(f'{purpose} 계산은 밀도 체크포인트 재시작을 지원하지 않습니다. 새 실행으로 시작하세요.')
    save_dft_config(config, directory / 'settings.json')
    scale = 1.0 if config.unit == 'angstrom' else BOHR_ANGSTROM
    initial = tuple(Atom(a.element, tuple(x * scale for x in a.position)) for a in config.atoms)
    # Intermediate steps need energy and forces only; keep their density cubes small.
    # Constraints belong to geomeTRIC; each geometry is a plain energy+gradient run.
    step_config = replace(config, atoms=initial, unit='angstrom', calculation='gradient', cube_points=12,
                          opt_constraints='')
    return directory, initial, step_config


def _geometric_molecule(geometric, atoms):
    molecule = geometric.molecule.Molecule()
    molecule.elem = [a.element for a in atoms]
    molecule.xyzs = [np.array([a.position for a in atoms])]
    return molecule


def optimize_molecule(config, directory, progress=print):
    try:
        import geometric
        import geometric.engine
        import geometric.molecule
        import geometric.optimize
        from geometric.errors import GeomOptNotConvergedError
    except ImportError as exc:
        raise RuntimeError('구조 최적화에는 geomeTRIC이 필요합니다. '
                           'wsl --exec .venv-dft/bin/python -m pip install -r requirements-optimize.txt') from exc
    directory, initial, step_config = _start(config, directory, '구조 최적화')
    criteria = OPT_CONVERGENCE[config.opt_convergence]
    molecule = _geometric_molecule(geometric, initial)
    engine, trajectory, state = _recording_engine(geometric.engine.Engine, molecule, step_config, directory,
                                                  progress, '구조 최적화')
    transition = config.calculation == 'transition_state'
    extra = {}
    constraints = config.opt_constraints.replace(';', '\n').strip()
    if constraints:
        from geometric.prepare import parse_constraints
        try:
            parse_constraints(molecule, constraints)
        except Exception as exc:  # geomeTRIC raises several types for malformed input.
            raise ValueError(f'geomeTRIC 제약 입력을 해석할 수 없습니다: {exc}') from exc
        (directory / 'constraints.txt').write_text(constraints + '\n', encoding='utf-8')
        extra['constraints'] = str(directory / 'constraints.txt')
    if transition:
        # P-RFO toward a first-order saddle; the starting Hessian comes from finite
        # differences of this engine's analytic gradients (calls recorded as kind=hessian).
        extra.update(transition=True, hessian='first')
    progress(f'geomeTRIC {geometric.__version__} {"전이상태 탐색" if transition else "구조 최적화"} · TRIC 내부좌표 · '
             f'기준 {config.opt_convergence} · 최대 {config.opt_max_steps}단계'
             + (' · 제약 있음' if constraints else ''))
    converged = True
    try:
        geometric.optimize.run_optimizer(
            customengine=engine, input=str(directory / 'geometric'), maxiter=config.opt_max_steps,
            logIni=str(Path(__file__).with_name('geometric_log.ini')),
            convergence_energy=criteria['energy'], convergence_grms=criteria['grms'],
            convergence_gmax=criteria['gmax'], convergence_drms=criteria['drms'], convergence_dmax=criteria['dmax'],
            **extra)
    except GeomOptNotConvergedError:
        converged = False
    if not trajectory:
        raise RuntimeError('geomeTRIC이 계산을 한 번도 요청하지 않았습니다. geometric.log를 확인하세요.')
    # The last evaluated geometry is the one whose energy/forces met (or failed) the criteria.
    final_atoms = state['atoms']
    progress('최종 구조를 요청한 출력 설정으로 다시 계산합니다.')
    final_config = replace(config, atoms=final_atoms, unit='angstrom', calculation='gradient', opt_constraints='')
    summary = calculate_molecule(final_config, directory, progress, state['density'])
    save_dft_config(final_config, directory / 'final_settings.json')
    save_dft_config(config, directory / 'settings.json')  # The request, not the final geometry.
    (directory / 'optimization.xyz').write_text(''.join(
        _xyz([Atom(a.element, tuple(p)) for a, p in zip(initial, row['coordinates_angstrom'])],
             f'step {row["step"]} E={row["energy_hartree"]:.12f} Ha') for row in trajectory), encoding='utf-8')
    (directory / 'optimized.xyz').write_text(
        _xyz(final_atoms, f'{config.formula} optimized={converged} E={summary["total_energy_hartree"]:.12f} Ha'),
        encoding='utf-8')
    gradient = np.asarray(summary['gradient']['gradient_hartree_per_bohr']) if summary['gradient'] else None
    # Under constraints the constrained coordinates keep a nonzero gradient by design.
    met = None if gradient is None or constraints else {
        'max_gradient': bool(np.abs(gradient).max() < criteria['gmax']),
        'rms_gradient': bool(np.sqrt(np.mean(gradient**2)) < criteria['grms'])}
    summary['optimization'] = {
        'optimizer': 'geomeTRIC', 'optimizer_version': geometric.__version__, 'coordinates': 'TRIC',
        'reference': 'L.-P. Wang, C. Song, J. Chem. Phys. 144, 214108 (2016)',
        'converged': converged, 'steps': len(trajectory), 'max_steps': config.opt_max_steps,
        'criteria_name': config.opt_convergence,
        'criteria': {'energy_change_hartree': criteria['energy'], 'rms_gradient_hartree_per_bohr': criteria['grms'],
                     'max_gradient_hartree_per_bohr': criteria['gmax'], 'rms_step_angstrom': criteria['drms'],
                     'max_step_angstrom': criteria['dmax']},
        'final_gradient_meets_force_criteria': met,
        'initial_energy_hartree': trajectory[0]['energy_hartree'],
        'final_energy_hartree': summary['total_energy_hartree'],
        # Same geometry as the last step, recomputed with the requested output settings.
        'final_recompute_energy_difference_hartree': summary['total_energy_hartree'] - trajectory[-1]['energy_hartree'],
        'energy_change_ev': (summary['total_energy_hartree'] - trajectory[0]['energy_hartree']) * HARTREE_TO_EV,
        'trajectory': trajectory, 'trajectory_file': 'optimization.xyz', 'optimized_xyz': 'optimized.xyz',
        'final_settings_file': 'final_settings.json', 'log_file': 'geometric.log',
        'constraints': constraints or None, 'transition_state_search': transition,
        'hessian_gradient_calculations': sum(row['kind'] == 'hessian' for row in trajectory),
        'completed_at_utc': datetime.now(timezone.utc).isoformat()}
    summary['calculation'] = config.calculation
    summary['config'] = config.to_dict()
    summary['final_geometry_config'] = final_config.to_dict()
    summary['assumptions'] = ['isolated finite molecule',
                              'first-order saddle (P-RFO) search near the input geometry' if transition else
                              'constrained local minimum search' if constraints else
                              'local minimum search from the input geometry; not a global or conformer search',
                              *summary['assumptions'][2:]]
    if transition:
        summary['warnings'].append('전이상태 후보입니다. 진동수 계산으로 허수 진동이 정확히 1개이고 그 모드가 '
                                   '원하는 반응 좌표인지 확인하세요(IRC는 미구현).')
    if config.grid_level < 3:
        summary['warnings'].insert(0, f'격자 수준 {config.grid_level}: 격자 응답 항이 빠진 힘의 오차가 수렴 기준과 '
                                      '비슷하거나 큽니다. 구조 최적화에는 grid_level 3 이상을 권장합니다.')
    if not converged:
        summary['warnings'].insert(0, f'구조 최적화 미수렴: {config.opt_max_steps}단계 안에 기준을 만족하지 못했습니다. '
                                      '최종 구조는 최소점이 아닙니다.')
    if not transition:
        summary['warnings'].append('최적화 구조는 입력 구조에서 출발한 국소 최소점 후보입니다. 진동수 분석(허수 진동 없음) 없이는 '
                               '최소점임을 확인하지 않았습니다.')
    scf_ok = summary['status'] == 'converged'
    summary['converged'] = scf_ok and converged
    summary['status'] = 'converged' if summary['converged'] else 'not_converged'
    temporary = directory / 'results.json.tmp'
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(directory / 'results.json')
    shutil.rmtree(directory / 'geometric.tmp', ignore_errors=True)
    progress('구조 최적화 수렴 · 결과 저장 완료' if summary['converged'] else '구조 최적화 미수렴 · 진단 결과 저장 완료')
    return summary


def scan_molecule(config, directory, progress=print):
    """Relaxed scan: one constrained geomeTRIC optimization per value, each started from the
    previous point's optimized geometry. Records live in scan/point_NNN."""
    kind, atoms, values = parse_scan(config.scan_spec, len(config.atoms))
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    save_dft_config(config, directory / 'settings.json')
    base = config.opt_constraints.replace(';', '\n').strip()
    rows, current, last = [], config, None
    for index, value in enumerate(values, 1):
        check_cancelled(directory)
        target = f'{kind} {" ".join(map(str, atoms))} {value:.10g}'
        point = replace(current, calculation='optimize', scan_spec='',
                        opt_constraints=(base + '\n' if base else '') + '$set\n' + target)
        progress(f'스캔 {index}/{len(values)} · {target}')
        last = optimize_molecule(point, directory / 'scan' / f'point_{index:03d}',
                                 lambda message: None if message.startswith(('OPT_JSON', 'SCF_JSON')) else progress(message))
        final = last['final_geometry_config']
        rows.append({'point': index, 'value': value, 'energy_hartree': last['total_energy_hartree'],
                     'converged': last['converged'], 'steps': last['optimization']['steps'],
                     'directory': f'scan/point_{index:03d}', 'coordinates_angstrom': [a['position'] for a in final['atoms']]})
        current = replace(current, atoms=tuple(Atom(**a) for a in final['atoms']), unit='angstrom')
    minimum = min(rows, key=lambda row: row['energy_hartree'])
    with (directory / 'scan.csv').open('w', encoding='utf-8') as stream:
        stream.write(f'point,{kind}_{"_".join(map(str, atoms))},energy_hartree,relative_kj_mol,converged\n')
        for row in rows:
            relative = (row['energy_hartree'] - minimum['energy_hartree']) * 2625.4996394799
            stream.write(f'{row["point"]},{row["value"]:.10g},{row["energy_hartree"]:.12f},{relative:.6f},{row["converged"]}\n')
    (directory / 'scan.xyz').write_text(''.join(
        _xyz([Atom(a.element, tuple(p)) for a, p in zip(config.atoms, row['coordinates_angstrom'])],
             f'point {row["point"]} {kind}={row["value"]:.6g} E={row["energy_hartree"]:.12f} Ha') for row in rows),
        encoding='utf-8')
    # The run folder shows the last point's density/orbitals, like any other run.
    for item in (directory / 'scan' / f'point_{len(rows):03d}').iterdir():
        if item.is_file() and item.name not in ('results.json', 'settings.json', 'final_settings.json'):
            shutil.copyfile(item, directory / item.name)
    summary = dict(last)
    summary['scan'] = {'coordinate': kind, 'atoms': list(atoms), 'units': 'angstrom' if kind == 'distance' else 'degree',
                       'points': rows, 'minimum_point': minimum['point'], 'table_file': 'scan.csv', 'xyz_file': 'scan.xyz',
                       'method': 'relaxed scan: constrained geomeTRIC optimization at each value'}
    summary['calculation'] = 'scan'
    summary['config'] = config.to_dict()
    summary['converged'] = all(row['converged'] for row in rows)
    summary['status'] = 'converged' if summary['converged'] else 'not_converged'
    summary['warnings'] = (['스캔의 일부 점이 수렴하지 않았습니다. scan.csv의 converged 열을 확인하세요.']
                           if not summary['converged'] else []) + [
        '최상위 결과(밀도·오비탈)는 마지막 스캔 점입니다. 이완 스캔 에너지는 반응 장벽의 근사이며 전이상태가 아닙니다.']
    temporary = directory / 'results.json.tmp'
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(directory / 'results.json')
    progress(f'스캔 완료 · 최저 에너지 점 {minimum["point"]} ({kind} = {minimum["value"]:.6g})')
    return summary
