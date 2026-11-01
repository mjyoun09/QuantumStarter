"""Fixed-geometry FFT mesh checks with independently retained DFT runs."""
import argparse
import csv
from dataclasses import replace
from datetime import datetime
import json
import math
from pathlib import Path
from tempfile import mkdtemp

import numpy as np

from dft_workflow import ROOT, CalculationError, run_calculation
from molecular import save_dft_config
from periodic import PeriodicConfig, load_calculation_config


def parse_meshes(text):
    """Parse '21 21 21;25 25 25;29 29 29'."""
    return tuple(tuple(int(x) for x in part.split()) for part in text.split(';'))


def validate_mesh_request(config, meshes, tolerances):
    if not isinstance(config, PeriodicConfig) or config.boundary != 'bulk' or config.integral_scheme != 'fftdf':
        raise ValueError('FFT mesh 수렴 확인은 bulk / fftdf 설정에서 지원합니다.')
    if len(meshes) < 3:
        raise ValueError('연속 두 구간 비교를 위해 FFT mesh를 3개 이상 지정하세요.')
    for index, mesh in enumerate(meshes):
        if len(mesh) != 3:
            raise ValueError('각 FFT mesh에 정수 3개가 필요합니다.')
        replace(config, fft_mesh=mesh, calculation='gradient')  # reuse config validation
        if index and (any(b < a for a, b in zip(meshes[index - 1], mesh)) or mesh == meshes[index - 1]):
            raise ValueError('FFT mesh는 성분별로 감소하지 않고 매 단계 하나 이상 증가해야 합니다.')
    if any(not math.isfinite(t) or t <= 0 for t in tolerances):
        raise ValueError('수렴 허용오차는 유한한 양수여야 합니다.')


def run_mesh_convergence(config, meshes, output_root=None, backend='auto', progress=print,
                         control=None, energy_tol_ev_per_atom=.001, force_tol_ev_a=.01,
                         stress_tol_gpa=.1):
    """Require the last two adjacent mesh comparisons to pass all three tolerances.

    This tests one geometry at fixed k points, basis and smearing, not completeness
    of the basis or Brillouin-zone sampling. F is compared when smearing is active.
    """
    meshes = tuple(tuple(mesh) for mesh in meshes)
    tolerances = (energy_tol_ev_per_atom, force_tol_ev_a, stress_tol_gpa)
    validate_mesh_request(config, meshes, tolerances)
    output_root = Path(output_root or ROOT / 'outputs' / 'fft_convergence').resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    directory = Path(mkdtemp(prefix=datetime.now().strftime('%Y%m%d_%H%M%S_%f_'), dir=output_root))
    save_dft_config(config, directory / 'settings.json')
    objective = 'free_energy_ev' if config.smearing_ev else 'total_energy_ev'
    rows, previous = [], None
    report = {'status': 'running', 'converged': False, 'config': config.to_dict(),
              'objective': objective, 'meshes': [list(m) for m in meshes], 'rows': rows,
              'tolerances': dict(zip(('energy_ev_per_atom', 'force_ev_per_angstrom', 'stress_gpa'), tolerances)),
              'rule': 'All three differences must pass for the last two consecutive mesh intervals.',
              'selected_mesh': None,
              'limitations': 'Fixed input geometry and scaled k mesh. A sampled numerical check, not proof of convergence; '
                             'repeat at the final relaxed geometry and check basis/k-point convergence separately.'}

    def save():
        (directory / 'mesh_convergence.json').write_text(
            json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
        fields = ['mesh', 'status', 'objective_ev', 'delta_energy_ev_per_atom',
                  'max_force_difference_ev_per_angstrom', 'max_stress_difference_gpa', 'passed', 'directory']
        with (directory / 'mesh_convergence.csv').open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)

    save()
    try:
        for mesh in meshes:
            if control is not None and control.requested:
                report['status'] = 'cancelled'
                break
            progress(f'FFT mesh 수렴 확인: {mesh}')
            report['current_mesh'] = list(mesh)
            save()
            cfg = replace(config, calculation='gradient', fft_mesh=mesh, band_path='')
            summary, run_dir = run_calculation(cfg, directory / 'runs', backend, progress, control=control)
            row = {'mesh': list(mesh), 'directory': str(Path(run_dir).relative_to(directory)),
                   'status': summary['status'], 'passed': False}
            rows.append(row)
            if not summary['converged'] or not summary.get('gradient'):
                report['status'] = 'cancelled' if summary['status'] == 'cancelled' else 'not_converged'
                save()
                break
            row['objective_ev'] = summary[objective]
            forces = np.asarray(summary['gradient']['forces_ev_per_angstrom'], dtype=float)
            stress = np.asarray(summary['gradient']['stress_gpa'], dtype=float)
            if not all(np.isfinite(x).all() for x in (forces, stress, np.asarray(row['objective_ev']))):
                raise ValueError('비유한 에너지/힘/응력은 수렴으로 판정할 수 없습니다.')
            if previous is not None:
                differences = (abs(row['objective_ev'] - previous[0]) / len(config.atoms),
                               float(np.max(np.linalg.norm(forces - previous[1], axis=1))),
                               float(np.max(np.abs(stress - previous[2]))))
                row.update(zip(('delta_energy_ev_per_atom', 'max_force_difference_ev_per_angstrom',
                                'max_stress_difference_gpa'), differences))
                row['passed'] = all(change <= tol for change, tol in zip(differences, tolerances))
            previous = (row['objective_ev'], forces, stress)
            save()
        else:
            report['converged'] = rows[-1]['passed'] and rows[-2]['passed']
            report['status'] = 'converged' if report['converged'] else 'not_converged'
            if report['converged']:
                report['selected_mesh'] = rows[-1]['mesh']
        if control is not None and control.requested:
            report.update(status='cancelled', converged=False, selected_mesh=None)
    except Exception as exc:
        report.update(status='error', error=str(exc), converged=False, selected_mesh=None)
        save()
        raise CalculationError(f'{exc}\nFFT 검사 기록: {directory}', directory) from exc
    save()
    return report, directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--meshes', default='21 21 21;25 25 25;29 29 29')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--backend', choices=['auto', 'local', 'wsl'], default='auto')
    parser.add_argument('--energy-tol', type=float, default=.001, help='eV/atom')
    parser.add_argument('--force-tol', type=float, default=.01, help='eV/Å')
    parser.add_argument('--stress-tol', type=float, default=.1, help='GPa')
    args = parser.parse_args()
    try:
        report, directory = run_mesh_convergence(load_calculation_config(args.config), parse_meshes(args.meshes),
                                                args.output, args.backend, energy_tol_ev_per_atom=args.energy_tol,
                                                force_tol_ev_a=args.force_tol, stress_tol_gpa=args.stress_tol)
        print(f"{report['status']} · 결과: {directory}")
        return 3 if report['status'] == 'cancelled' else (0 if report['converged'] else 2)
    except (ValueError, RuntimeError, OSError) as exc:
        print(f'실행 오류: {exc}')
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
