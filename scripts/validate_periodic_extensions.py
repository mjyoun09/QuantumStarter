"""Stage 10 probes: 2D bands and bulk-only FFTDF derivatives (not a product API)."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import time
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def slab(root):
    from ase.build import graphene
    from periodic import PeriodicConfig
    from engines import calculate
    cfg = PeriodicConfig.from_ase(graphene(a=2.46, vacuum=8), basis='gth-szv-molopt-sr',
                                  kmesh=(2, 2, 1), grid_level=2, precision=1e-7,
                                  conv_tol=1e-10, smearing_ev=.05, band_path='GMKG',
                                  band_points=16, threads=1)
    result = calculate(cfg, root / 'graphene', print)
    bands = result['band_structure']
    checks = {'converged': result['converged'], 'pbc': bands['pbc'],
              'max_abs_kz': float(np.max(np.abs(np.asarray(bands['kpoints_scaled'])[:, 2]))),
              'mesh_consistency_max_ev': bands['mesh_consistency_max_ev'],
              'states_error': result['dos']['integrated_states'] - result['dos']['expected_states'],
              'electron_error': result['electron_trace'] - result['explicit_electrons']}
    return checks


def fftdf(root):
    from pyscf import lib
    from pyscf.pbc import dft
    from periodic import load_calculation_config
    from engines.periodic import build_cell
    lib.num_threads(1)
    cfg = replace(load_calculation_config('examples/silicon_bulk.json'), basis='gth-szv-molopt-sr',
                  kmesh=(1, 1, 1), precision=1e-8, conv_tol=1e-11, threads=1)
    cfg = replace(cfg, atoms=(cfg.atoms[0], replace(cfg.atoms[1], position=(1.4375, 1.3575, 1.3575))))
    base, _ = build_cell(cfg)
    mesh = (25, 25, 25)
    rows = []

    def solve(lattice, positions, label):
        cell = base.copy()
        cell.a, cell.atom, cell.unit = lattice, [('Si', p) for p in positions], 'bohr'
        cell.mesh = np.asarray(mesh)
        cell.build(False, False)
        mf = dft.KRKS(cell, cell.make_kpts(cfg.kmesh))  # default FFTDF, no density_fit()
        mf.xc, mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = 'pbe', 1e-11, 1e-7, 100
        mf.small_rho_cutoff = 0
        mf.grids = dft.gen_grid.UniformGrids(cell)
        mf.grids.mesh = np.asarray(mesh)
        energy = float(mf.kernel())
        row = {'label': label, 'energy_hartree': energy, 'converged': bool(mf.converged),
               'mesh': list(mesh), 'lattice_bohr': lattice.tolist(), 'positions_bohr': positions.tolist()}
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if k not in ('lattice_bohr', 'positions_bohr')}), flush=True)
        if not mf.converged:
            raise RuntimeError(f'{label}: SCF did not converge')
        return mf

    lattice, positions = base.lattice_vectors(), base.atom_coords()
    mf = solve(lattice, positions, 'reference')
    grad = mf.nuc_grad_method()
    gradient = grad.kernel()
    stress = grad.get_stress()
    h = .001
    delta = np.zeros_like(positions)
    delta[1, 0] = h
    finite_gradient = (solve(lattice, positions + delta, 'atom_plus').e_tot -
                       solve(lattice, positions - delta, 'atom_minus').e_tot) / (2 * h)
    strain_step = 1e-4
    finite_stress = {}
    for x, y in ((0, 0), (1, 2)):
        delta = np.zeros((3, 3))
        delta[x, y] = strain_step
        plus, minus = np.eye(3) + delta, np.eye(3) - delta
        ep = solve(lattice @ plus, positions @ plus, f'strain_{x}{y}_plus').e_tot
        em = solve(lattice @ minus, positions @ minus, f'strain_{x}{y}_minus').e_tot
        finite_stress[f'{x}{y}'] = (ep - em) / (2 * strain_step * base.vol)
    report = {'scope': 'bulk KRKS/PBE/GTH, gamma, uniform fixed FFT mesh; no smearing',
              'pyscf_stress_driver': 'mixed analytic and internal finite differences',
              'gradient_hartree_per_bohr': gradient.tolist(),
              'finite_gradient_atom1_x': finite_gradient,
              'gradient_error': float(gradient[1, 0] - finite_gradient),
              'stress_hartree_per_bohr3': stress.tolist(), 'finite_stress': finite_stress,
              'stress_error': {xy: float(stress[int(xy[0]), int(xy[1])] - value)
                               for xy, value in finite_stress.items()}, 'runs': rows}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['slab', 'fftdf'])
    parser.add_argument('--label', required=True)
    args = parser.parse_args()
    root = Path('outputs/periodic_extension_validation') / args.label / args.mode
    root.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    result = (slab if args.mode == 'slab' else fftdf)(root)
    result['seconds'] = time.perf_counter() - start
    (root / 'report.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k != 'runs'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
