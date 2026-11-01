"""Validate bulk cell relaxation and energies at nearby volumes (numerical test only)."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from engines import calculate
from periodic import load_calculation_config, PeriodicConfig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', required=True)
    parser.add_argument('--mode', choices=['full', 'isotropic'], default='isotropic')
    parser.add_argument('--steps', type=int, default=40)
    parser.add_argument('--shear', type=float, default=0., help='Initial symmetric xy strain')
    parser.add_argument('--smearing', type=float, default=.1, help='Fermi width in eV; zero tests E')
    parser.add_argument('--check-volumes', action='store_true')
    args = parser.parse_args()
    root = Path('outputs/cell_validation') / args.label
    root.mkdir(parents=True, exist_ok=False)
    cfg = replace(load_calculation_config('examples/silicon_cell.json'), cell_relaxation=args.mode,
                  cell_max_steps=args.steps, smearing_ev=args.smearing)
    if args.shear:
        deformation = np.eye(3)
        deformation[0, 1] = deformation[1, 0] = args.shear
        cfg = replace(cfg, lattice=tuple(map(tuple, np.asarray(cfg.lattice) @ deformation)),
                      atoms=tuple(replace(a, position=tuple(np.asarray(a.position) @ deformation)) for a in cfg.atoms))
    result = calculate(cfg, root / 'relax')
    report = {'cell_optimization': result['cell_optimization'], 'gradient': result['gradient'],
              'elapsed_seconds': result['elapsed_seconds'], 'config': cfg.to_dict(),
              'final_geometry_config': result['final_geometry_config']}
    if args.check_volumes and result['converged']:
        final = PeriodicConfig.from_dict(result['final_geometry_config'])
        objective = 'free_energy_ev' if final.smearing_ev else 'total_energy_ev'
        changes = []
        for scale in (.995, 1.005):
            displaced = replace(final, calculation='single_point',
                                lattice=tuple(map(tuple, np.asarray(final.lattice) * scale)),
                                atoms=tuple(replace(a, position=tuple(np.asarray(a.position) * scale)) for a in final.atoms))
            trial = calculate(displaced, root / f'volume_{scale}', lambda _: None)
            if not trial['converged']:
                raise RuntimeError('Volume comparison SCF did not converge')
            changes.append({'linear_scale': scale, 'delta_objective_ev': trial[objective] - result[objective]})
        report['nearby_volume_changes'] = changes
    (root / 'validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)
    if not result['converged']:
        raise RuntimeError('Cell relaxation did not converge; retained as a nonconverged run')
    if any(row['delta_objective_ev'] <= 0 for row in report.get('nearby_volume_changes', [])):
        raise RuntimeError('A nearby sampled volume has a lower objective')


if __name__ == '__main__':
    main()
