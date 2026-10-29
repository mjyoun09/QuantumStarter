"""Persist native/reference runs under identical conditions; no reference energies hardcoded."""
from dataclasses import replace
from datetime import datetime
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular import DFTConfig, parse_coordinates, save_dft_config
from dft_workflow import run_calculation


def main():
    output = ROOT / 'outputs' / 'native_validation' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    rows = []
    cases = [('H', 'H 0 0 0', 1, 'lda,vwn'), ('He', 'He 0 0 0', 0, 'lda,vwn'), ('H2', 'H 0 0 0; H 0 0 .74', 0, 'lda,vwn'),
             ('H2O', 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587', 0, 'lda,vwn'),
             ('OH', 'O 0 0 0; H 0 0 .97', 1, 'pbe')]
    for name, atoms, spin, xc in cases:
        config = DFTConfig(atoms=parse_coordinates(atoms), spin=spin, xc=xc, engine='native',
                           basis='sto-3g', conv_tol=1e-10, threads=2, grid_level=2, cube_points=16)
        native, npath = run_calculation(config, output / name, 'local')
        external, epath = run_calculation(replace(config, engine='pyscf'), output / name, 'local')
        difference = native['total_energy_hartree'] - external['total_energy_hartree']
        assert native['converged'] and external['converged']
        aligned, apath = run_calculation(replace(config, engine='pyscf'), output / name, 'local', restart_from=npath)
        aligned_difference = native['total_energy_hartree'] - aligned['total_energy_hartree']
        assert aligned['converged'] and abs(aligned_difference) < 2e-7
        if spin == 0:
            assert abs(difference) < 2e-7
        rows.append({'case': name, 'native_hartree': native['total_energy_hartree'],
                     'pyscf_hartree': external['total_energy_hartree'], 'difference_hartree': difference,
                     'pyscf_from_native_density_hartree': aligned['total_energy_hartree'],
                     'aligned_difference_hartree': aligned_difference, 'aligned_directory': str(apath.relative_to(ROOT)),
                     'native_directory': str(npath.relative_to(ROOT)), 'pyscf_directory': str(epath.relative_to(ROOT)),
                     'native_final_residual': native['scf_criteria']['final_residual']})
        (output / 'comparison.json').write_text(json.dumps(rows, indent=2) + '\n', encoding='utf-8')
        print(rows[-1], flush=True)
    print('REPORT', output / 'comparison.json')


if __name__ == '__main__':
    main()
