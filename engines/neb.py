"""Climbing-image nudged elastic band between two fixed endpoints (geomeTRIC chain optimizer).

geomeTRIC 1.1 optimizes the chain (NEB forces, climbing image, endpoints locked);
every image energy and gradient comes from this project's engines. geomeTRIC does
not interpolate, so the initial path is a linear Cartesian interpolation between
the reactant (config.atoms) and the product (config.neb_product), with the atom
order and orientation exactly as given (no automatic alignment or reordering).
"""
from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import numpy as np

from molecular import Atom, parse_coordinates, save_dft_config
from scf_control import check_cancelled
from .optimizer import _xyz
from .pyscf_engine import BOHR_ANGSTROM, calculate_molecule

HARTREE_TO_KJ_MOL = 2625.4996394799


def neb_molecule(config, directory, progress=print):
    try:
        import geometric
        import geometric.engine
        import geometric.molecule
        from geometric.neb import ElasticBand, OptimizeChain
        from geometric.params import NEBParams
    except ImportError as exc:
        raise RuntimeError('NEB에는 geomeTRIC이 필요합니다. '
                           'wsl --exec .venv-dft/bin/python -m pip install -r requirements-optimize.txt') from exc
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    save_dft_config(config, directory / 'settings.json')
    scale = 1.0 if config.unit == 'angstrom' else BOHR_ANGSTROM
    reactant = np.array([a.position for a in config.atoms]) * scale
    product = np.array([a.position for a in parse_coordinates(config.neb_product)]) * scale
    elements = [a.element for a in config.atoms]
    fractions = np.linspace(0, 1, config.neb_images)
    path = [reactant + f * (product - reactant) for f in fractions]
    base = replace(config, calculation='gradient', unit='angstrom', cube_points=12, neb_product='', opt_constraints='')
    calls, cycles, densities = [], [], {}

    def quiet(message):
        check_cancelled(directory)
        if message.startswith(('계산 오류', 'SCF 미수렴')):
            progress(message)

    class Engine(geometric.engine.Engine):
        def calc_new(self, coords, dirname):
            check_cancelled(directory)
            index = len(calls) + 1
            positions = np.asarray(coords, dtype=float).reshape(-1, 3) * BOHR_ANGSTROM
            atoms = tuple(Atom(e, tuple(p)) for e, p in zip(elements, positions))
            image = Path(str(dirname)).name  # geomeTRIC names image folders struct_NN
            target = directory / 'neb' / f'calc_{index:04d}'
            summary = calculate_molecule(replace(base, atoms=atoms), target, quiet, densities.get(image))
            if not summary['converged'] or summary['gradient'] is None:
                raise RuntimeError(f'NEB 이미지 계산 {target}의 SCF가 수렴하지 않아 중단합니다.')
            with np.load(target / 'electronic_structure.npz', allow_pickle=False) as data:
                densities[image] = np.asarray(data['density_matrix'])
            calls.append({'call': index, 'image': image, 'energy_hartree': summary['total_energy_hartree'],
                          'directory': f'neb/calc_{index:04d}'})
            return {'energy': summary['total_energy_hartree'],
                    'gradient': np.asarray(summary['gradient']['gradient_hartree_per_bohr']).ravel()}

    molecule = geometric.molecule.Molecule()
    molecule.elem = elements
    molecule.xyzs = path
    params = NEBParams(neb=True, images=config.neb_images, maxg=config.neb_maxg_ev_per_angstrom,
                       avgg=config.neb_maxg_ev_per_angstrom / 2, neb_maxcyc=config.opt_max_steps,
                       prefix=str(directory / 'geometric_neb'))
    tmpdir = directory / 'geometric_neb.tmp'
    tmpdir.mkdir(exist_ok=True)

    def record(chain):
        energies = [float(s.energy) for s in chain.Structures]
        cycles.append({'cycle': len(cycles) + 1, 'energies_hartree': energies,
                       'max_gradient_ev_per_angstrom': float(chain.maxg) if hasattr(chain, 'maxg') else None,
                       'avg_gradient_ev_per_angstrom': float(chain.avgg) if hasattr(chain, 'avgg') else None})
        barrier = (max(energies[1:-1]) - energies[0]) * HARTREE_TO_KJ_MOL if len(energies) > 2 else 0.0
        progress(f'NEB {len(cycles):3d}회 · 최고 이미지 장벽 {barrier:.3f} kJ/mol · '
                 + (f'최대 힘 {cycles[-1]["max_gradient_ev_per_angstrom"]:.4f} eV/Å'
                    if cycles[-1]['max_gradient_ev_per_angstrom'] is not None else ''))

    progress(f'geomeTRIC {geometric.__version__} CI-NEB · 이미지 {config.neb_images}개(양 끝 고정) · '
             f'수렴 최대/평균 {params.maxg}/{params.avgg} eV/Å · 최대 {config.opt_max_steps}회')
    import logging.config
    logging.config.fileConfig(str(Path(__file__).with_name('geometric_log.ini')),
                              defaults={'logfilename': str(directory / 'geometric_neb.log')},
                              disable_existing_loggers=False)
    engine = Engine(molecule)
    chain = ElasticBand(molecule, engine=engine, tmpdir=str(tmpdir), params=params, plain=0)
    chain, iterations = OptimizeChain(chain, engine, params, save_callback=record)
    maxg, avgg = float(chain.maxg), float(chain.avgg)
    converged = maxg < params.maxg and avgg < params.avgg
    images = [{'image': i, 'energy_hartree': float(s.energy),
               'coordinates_angstrom': (np.asarray(s.cartesians).reshape(-1, 3) * BOHR_ANGSTROM).tolist()}
              for i, s in enumerate(chain.Structures)]
    reference = images[0]['energy_hartree']
    for row in images:
        row['relative_kj_mol'] = (row['energy_hartree'] - reference) * HARTREE_TO_KJ_MOL
    top = max(images[1:-1], key=lambda row: row['energy_hartree'])
    (directory / 'neb_path.xyz').write_text(''.join(
        _xyz([Atom(e, tuple(p)) for e, p in zip(elements, row['coordinates_angstrom'])],
             f'image {row["image"]} E={row["energy_hartree"]:.12f} Ha rel={row["relative_kj_mol"]:.4f} kJ/mol')
        for row in images), encoding='utf-8')
    with (directory / 'neb.csv').open('w', encoding='utf-8') as stream:
        stream.write('image,energy_hartree,relative_kj_mol\n')
        stream.writelines(f'{r["image"]},{r["energy_hartree"]:.12f},{r["relative_kj_mol"]:.6f}\n' for r in images)
    progress(f'최고 에너지 이미지 {top["image"]}를 요청한 출력 설정으로 다시 계산합니다.')
    top_atoms = tuple(Atom(e, tuple(p)) for e, p in zip(elements, top['coordinates_angstrom']))
    final_config = replace(base, atoms=top_atoms, cube_points=config.cube_points)
    summary = calculate_molecule(final_config, directory, progress)
    save_dft_config(config, directory / 'settings.json')
    save_dft_config(final_config, directory / 'final_settings.json')
    summary['neb'] = {
        'method': 'climbing-image NEB (geomeTRIC chain optimizer), endpoints fixed',
        'optimizer_version': geometric.__version__, 'images': len(images), 'cycles': len(cycles),
        'iterations_reported': int(iterations) if isinstance(iterations, (int, np.integer)) else None,
        'converged': converged, 'max_gradient_ev_per_angstrom': maxg, 'avg_gradient_ev_per_angstrom': avgg,
        'criteria_ev_per_angstrom': {'max': params.maxg, 'average': params.avgg},
        'spring_constant_kcal_mol_angstrom2': params.nebk, 'initial_path': 'linear Cartesian interpolation, no alignment',
        'image_data': images, 'highest_image': top['image'],
        'barrier_from_reactant_kj_mol': top['relative_kj_mol'],
        'reaction_energy_kj_mol': images[-1]['relative_kj_mol'],
        'gradient_calculations': len(calls), 'cycle_history': cycles,
        'path_file': 'neb_path.xyz', 'table_file': 'neb.csv', 'log_file': 'geometric_neb.log',
        'completed_at_utc': datetime.now(timezone.utc).isoformat()}
    summary['calculation'] = 'neb'
    summary['config'] = config.to_dict()
    summary['final_geometry_config'] = final_config.to_dict()
    summary['assumptions'] = ['isolated finite molecule', 'minimum energy path between the given fixed endpoints',
                              *summary['assumptions'][2:]]
    summary['warnings'] = ([] if converged else [
        f'NEB 미수렴: 최대/평균 힘 {maxg:.4f}/{avgg:.4f} eV/Å가 기준({params.maxg}/{params.avgg})을 넘습니다.']) + [
        '최고 에너지 이미지는 전이상태 추정입니다. transition_state 계산으로 정밀화하고 진동수로 허수 진동 1개를 확인하세요.',
        '양 끝 구조는 고정되며 최적화되지 않습니다. 각 끝이 최소점인지 먼저 확인하세요.'] + summary['warnings']
    summary['converged'] = converged and summary['status'] == 'converged'
    summary['status'] = 'converged' if summary['converged'] else 'not_converged'
    temporary = directory / 'results.json.tmp'
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(directory / 'results.json')
    shutil.rmtree(tmpdir, ignore_errors=True)
    progress(f'NEB {"수렴" if converged else "미수렴"} · 장벽 {top["relative_kj_mol"]:.3f} kJ/mol '
             f'({top["relative_kj_mol"] / 96.4853321233 :.4f} eV) · 이미지 {top["image"]}')
    return summary
