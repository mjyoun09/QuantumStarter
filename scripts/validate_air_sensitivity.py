"""Propagate a two-basis, two-grid PBE/RRHO comparison into the N2/O2 mixture."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from engines.gas_thermo import GASES, write_csv
from gas_mixture import DEFAULT_MIX, mixture_report
from molecular import DFTConfig, parse_coordinates
from molecular_convergence import DEFAULT_BASES, DEFAULT_GRIDS, run_comparison


def run(out, backend):
    out.mkdir(parents=True, exist_ok=False)
    coefficients = {(b, g): [] for b in DEFAULT_BASES for g in DEFAULT_GRIDS}
    for species in ('N2', 'O2'):
        coords, spin, *_ = GASES[species]
        cfg = DFTConfig(calculation='gas_thermo', atoms=parse_coordinates(coords), spin=spin,
                        engine='native', xc='pbe', cube_points=12, threads=1, opt_convergence='tight',
                        conv_tol=1e-10, native_density_tol=1e-8, native_residual_tol=1e-8)
        report, directory = run_comparison(cfg, output_root=out/species, backend=backend)
        if report['status'] != 'completed':
            raise RuntimeError(f'{species}: {directory}/comparison.json')
        for row in report['runs']:
            coefficients[row['basis'], row['grid']].append(directory/row['directory']/'nasa7.yaml')
    table = []
    for (basis, grid), paths in coefficients.items():
        for row in mixture_report(DEFAULT_MIX, paths)['rows']:
            table.append(dict(basis=basis, grid=grid, **row))
    write_csv(out/'air_sensitivity.csv', table)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--backend', choices=('auto','local','wsl'), default='auto')
    args = parser.parse_args()
    run(args.out, args.backend)
