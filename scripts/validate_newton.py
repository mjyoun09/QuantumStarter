"""Reproducible native Newton probes; run with the WSL .venv-dft Python."""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular import DFTConfig, parse_coordinates


def main():
    from pyscf import gto, lib
    from engines.native_scf import NativeKS
    from engines.pyscf_engine import configure_grids

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['water', 'Fe', 'Cr', 'Co'], default='Cr')
    parser.add_argument('--max-cycle', type=int, default=100)
    parser.add_argument('--label', default='probe')
    parser.add_argument('--benchmark', action='store_true', help='Compare analytic and FD H.v at the same orbitals')
    args = parser.parse_args()
    root = Path('outputs/newton_validation') / args.label / args.case
    root.mkdir(parents=True, exist_ok=False)
    atoms, spin = ('O 0 0 0; H 0 -.757 .587; H 0 .757 .587', 0) if args.case == 'water' else (
        f'{args.case} 0 0 0', {'Fe': 4, 'Cr': 6, 'Co': 3}[args.case])
    cfg = DFTConfig(atoms=parse_coordinates(atoms), spin=spin, xc='pbe', basis='sto-3g',
                    engine='native', scf_solver='newton', grid_level=3, threads=1,
                    max_cycle=args.max_cycle)
    lib.num_threads(cfg.threads)
    mol = gto.M(atom=atoms, spin=spin, basis=cfg.basis, verbose=0)
    mf = NativeKS(mol, cfg, root)
    configure_grids(mf.grids, cfg)
    rows = []

    def callback(env):
        row = {key: env[key] for key in ('cycle', 'e_tot', 'norm_gorb', 'norm_ddm', 'delta_energy')}
        rows.append(row)
        if row['cycle'] % 10 == 0:
            print(json.dumps(row), flush=True)

    mf.callback = callback
    start = time.perf_counter()
    energy = mf.kernel()
    seconds = time.perf_counter() - start
    report = {'case': args.case, 'atoms': atoms, 'spin': spin, 'xc': cfg.xc, 'basis': cfg.basis,
              'grid_level': cfg.grid_level, 'threads': cfg.threads, 'max_cycle': cfg.max_cycle,
              'criteria': {'energy': cfg.conv_tol, 'density': cfg.native_density_tol,
                           'residual': cfg.native_residual_tol},
              'energy_hartree': energy, 'converged': mf.converged, 'seconds': seconds,
              'cycles': len(rows), 'final_residual': mf.final_residual,
              'aufbau_violation_hartree': mf.aufbau_violation_hartree,
              'curvature_checks': mf.newton_escapes, 'history': rows}
    if args.benchmark:
        from engines.native_stability import OrbitalRotations
        rotations = OrbitalRotations(mf)
        vector = np.random.default_rng(47).normal(size=rotations.dimension)
        vector /= np.linalg.norm(vector)
        prepared = time.perf_counter()
        analytic = rotations.hessian_vector(vector)
        setup_seconds = time.perf_counter() - prepared
        finite = rotations.hessian_vector_fd(vector)
        timings = {}
        for label, operation in [('analytic', rotations.hessian_vector), ('finite_difference', rotations.hessian_vector_fd)]:
            times = []
            for _ in range(7):
                started = time.perf_counter()
                operation(vector)
                times.append(time.perf_counter() - started)
            timings[label] = float(np.median(times))
        report['response_benchmark'] = {'repeats': 7, 'threads': cfg.threads,
            'first_analytic_call_seconds': setup_seconds, 'median_seconds': timings,
            'speedup': timings['finite_difference'] / timings['analytic'],
            'max_hessian_vector_difference': float(np.max(np.abs(analytic - finite)))}
    (root / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'history'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
