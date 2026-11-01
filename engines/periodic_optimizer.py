"""Bulk FFTDF cell/atom relaxation using ASE BFGS and FrechetCellFilter.

SCF, forces and stress are PySCF. ASE supplies the optimizer and logarithmic
cell coordinates. Fixed integer FFT mesh and scaled k mesh throughout; minimize
E, or the Mermin free energy F with smearing, at zero external pressure.
"""
import csv
from dataclasses import replace
import json
from pathlib import Path
import shutil
import time
import numpy as np
from scipy.constants import physical_constants

from molecular import save_dft_config
from periodic import PeriodicAtom
from scf_control import check_cancelled
from .periodic import build_cell, calculate_periodic

EV_PER_ANGSTROM3_GPA = physical_constants['electron volt'][0] / 1e-30 / 1e9


def stress_residual(stress_gpa, mode):
    stress = np.asarray(stress_gpa)
    symmetric = (stress + stress.T) / 2
    return float(abs(np.trace(symmetric) / 3) if mode == 'isotropic' else np.abs(symmetric).max())


def optimize_periodic(config, directory, progress=print):
    import ase
    from ase import Atoms
    from ase.calculators.calculator import Calculator, all_changes
    from ase.filters import FrechetCellFilter
    from ase.io import write
    from ase.optimize import BFGS

    if config.integral_scheme != 'fftdf' or config.boundary != 'bulk':
        raise ValueError('셀 최적화는 FFTDF bulk에서만 지원합니다.')
    started = time.perf_counter()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    save_dft_config(config, directory / 'settings.json')
    check_cancelled(directory)
    # Resolve an automatic mesh once. Rechoosing it while deforming the cell
    # would change the discretized objective discontinuously.
    mesh = config.fft_mesh or tuple(int(n) for n in build_cell(config)[0].mesh)
    base = replace(config, calculation='gradient', fft_mesh=mesh, band_path='')
    objective_key = 'free_energy_ev' if config.smearing_ev else 'total_energy_ev'
    atoms = Atoms([a.element for a in config.atoms], positions=[a.position for a in config.atoms],
                  cell=config.lattice, pbc=True)
    records = []
    state = {}

    def geometry_config(atoms):
        return replace(base, lattice=tuple(map(tuple, atoms.cell.array)),
                       atoms=tuple(PeriodicAtom(symbol, tuple(x)) for symbol, x in zip(atoms.symbols, atoms.positions)))

    class PeriodicCalculator(Calculator):
        implemented_properties = ['energy', 'free_energy', 'forces', 'stress']

        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            check_cancelled(directory)
            cfg = geometry_config(self.atoms)
            step_dir = directory / 'steps' / f'step_{len(records):04d}'
            summary = calculate_periodic(cfg, step_dir, progress, cancel_directory=directory)
            if not summary['converged'] or summary['gradient'] is None:
                raise RuntimeError(f'셀 최적화 중 SCF 미수렴: {step_dir}. 힘/응력을 최적화에 전달하지 않았습니다.')
            gradient = summary['gradient']
            stress = np.asarray(gradient['stress_gpa'])
            self.results = dict(energy=summary[objective_key], free_energy=summary[objective_key],
                                forces=np.asarray(gradient['forces_ev_per_angstrom']),
                                stress=(stress + stress.T) / (2 * EV_PER_ANGSTROM3_GPA))
            previous = records[-1]['objective_ev'] if records else None
            row = {'step': len(records), 'objective_ev': summary[objective_key],
                   'energy_change_ev': None if previous is None else summary[objective_key] - previous,
                   'volume_angstrom3': float(self.atoms.get_volume()),
                   'max_force_ev_per_angstrom': gradient['max_atomic_force_ev_per_angstrom'],
                   'stress_residual_gpa': stress_residual(stress, config.cell_relaxation),
                   'stress_gpa': stress.tolist(), 'lattice_angstrom': self.atoms.cell.array.tolist(),
                   'coordinates_angstrom': self.atoms.positions.tolist(),
                   'directory': str(step_dir.relative_to(directory)), 'scf_converged': True}
            records.append(row)
            state.update(config=cfg, summary=summary)
            (directory / 'cell_history.json').write_text(json.dumps(records, indent=2) + '\n', encoding='utf-8')
            write(directory / 'cell_optimization.extxyz', self.atoms, format='extxyz', append=len(records) > 1,
                  write_results=False)
            progress(f'셀 최적화 {row["step"]}: V={row["volume_angstrom3"]:.6f} Å³ · '
                     f'최대 힘 {row["max_force_ev_per_angstrom"]:.3e} eV/Å · 응력 잔차 {row["stress_residual_gpa"]:.3e} GPa')

    atoms.calc = PeriodicCalculator()
    cell_filter = FrechetCellFilter(atoms, hydrostatic_strain=config.cell_relaxation == 'isotropic')
    optimizer = BFGS(cell_filter, logfile=str(directory / 'cell_optimizer.log'), maxstep=.1)

    def criteria(row):
        return {'force': row['max_force_ev_per_angstrom'] <= config.cell_force_tol,
                'stress': row['stress_residual_gpa'] <= config.cell_stress_tol_gpa,
                'energy_change': row['energy_change_ev'] is None or abs(row['energy_change_ev']) <= config.cell_energy_tol_ev}

    # ASE's combined pseudo-force norm mixes atomic/cell coordinates. Use our
    # physical force/stress/energy criteria instead of accepting that norm alone.
    with optimizer:
        for _ in optimizer.irun(fmax=0., steps=config.cell_max_steps):
            check_cancelled(directory)
            if all(criteria(records[-1]).values()):
                break
    final_config = replace(state['config'], band_path=config.band_path)
    # Never publish a converged single-point result at the optimization root:
    # a crash/cancel during final verification must not masquerade as convergence.
    final_directory = directory / 'final'
    summary = calculate_periodic(final_config, final_directory, progress, cancel_directory=directory)
    check_cancelled(directory)
    final_gradient = summary['gradient']
    final_met = {'scf': bool(summary['converged']), 'force': False, 'stress': False, 'energy_change': False}
    if final_gradient is not None:
        final_met.update(force=final_gradient['max_atomic_force_ev_per_angstrom'] <= config.cell_force_tol,
                         stress=stress_residual(final_gradient['stress_gpa'], config.cell_relaxation) <= config.cell_stress_tol_gpa,
                         energy_change=criteria(records[-1])['energy_change'] and
                         abs(summary[objective_key] - records[-1]['objective_ev']) <= config.cell_energy_tol_ev)
    converged = all(final_met.values())
    for artifact in final_directory.iterdir():
        if artifact.is_file() and artifact.name not in ('settings.json', 'results.json'):
            shutil.copy2(artifact, directory / artifact.name)
    save_dft_config(config, directory / 'settings.json')
    save_dft_config(final_config, directory / 'final_settings.json')
    write(directory / 'optimized.extxyz', atoms, format='extxyz', write_results=False)
    fields = ['step', 'objective_ev', 'energy_change_ev', 'volume_angstrom3',
              'max_force_ev_per_angstrom', 'stress_residual_gpa', 'directory']
    with (directory / 'cell_optimization.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    summary.update(config=config.to_dict(), final_geometry_config=final_config.to_dict(),
                   calculation='cell_optimize', converged=converged,
                   status='converged' if converged else 'not_converged', elapsed_seconds=time.perf_counter() - started)
    summary['cell_optimization'] = {
        'optimizer': 'ASE BFGS / FrechetCellFilter', 'optimizer_version': ase.__version__,
        'mode': config.cell_relaxation, 'external_pressure_gpa': 0.,
        'objective': 'Mermin free energy F = E - TS' if config.smearing_ev else 'total energy E',
        'converged': converged, 'steps': optimizer.nsteps, 'max_steps': config.cell_max_steps,
        'criteria': {'force_ev_per_angstrom': config.cell_force_tol, 'stress_gpa': config.cell_stress_tol_gpa,
                     'energy_change_ev_per_cell': config.cell_energy_tol_ev}, 'final_criteria_met': final_met,
        'stress_residual_definition': 'abs(trace(sigma)/3)' if config.cell_relaxation == 'isotropic' else 'max(abs(sym(sigma)))',
        'fft_mesh_policy': 'fixed integer mesh from initial cell', 'fft_mesh': list(mesh),
        'kpoint_policy': 'fixed scaled k mesh; Cartesian k points change with cell',
        'history_file': 'cell_history.json', 'trajectory_file': 'cell_optimization.extxyz',
        'final_structure_file': 'optimized.extxyz', 'final_settings_file': 'final_settings.json',
        'final_single_point_directory': 'final',
        'initial_objective_ev': records[0]['objective_ev'], 'final_objective_ev': summary[objective_key],
        'final_recompute_difference_ev': summary[objective_key] - records[-1]['objective_ev']}
    summary['warnings'].append('셀 최적화는 지정된 FFT mesh와 k점에서의 국소 정지점 후보입니다. '
                               '최종 구조에서도 더 촘촘한 FFT mesh·k점 검증 및 동역학적 안정성 확인이 필요합니다.')
    if config.cell_relaxation == 'isotropic':
        summary['warnings'].append('등방 셀 최적화: 셀 모양을 고정하고 체적만 조절하므로 비등방 응력은 남을 수 있습니다.')
    if not converged:
        summary['warnings'].insert(0, '셀 최적화 미수렴: 최종 SCF와 원자 힘·응력·에너지 변화 기준을 확인하세요.')
    (directory / 'results.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    return summary
