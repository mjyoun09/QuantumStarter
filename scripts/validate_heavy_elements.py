"""Boundary-path checks; these are not ground-state or experimental benchmarks.

Run with the PySCF environment from the project root. Each execution gets new
run directories, retaining input, logs, electronic structure and provenance.
"""
from datetime import datetime
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from dft_workflow import run_calculation
from molecular import Atom, DFTConfig


def main():
    output = ROOT / 'outputs' / 'all_elements_validation' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    report = []
    for symbol, charge, basis, grid in [('Cs', 1, 'def2-svp', 2), ('Ce', 4, 'dyall-v2z', 2),
                                       ('U', 6, 'dyall-v2z', 2), ('Og', 0, 'dyall-v2z', 3)]:
        config = DFTConfig(atoms=(Atom(symbol, (0, 0, 0)),), charge=charge, basis=basis,
                           spin=0, scf_solver='newton', threads=4, memory_mb=4000,
                           grid_level=grid, cube_points=16, conv_tol=1e-8)
        result, directory = run_calculation(config, output / symbol, backend='local', progress=print)
        data = np.load(directory / 'electronic_structure.npz', allow_pickle=False)
        coeff, overlap = data['mo_coeff'], data['overlap']
        assert result['converged'], (symbol, 'not converged')
        assert abs(result['electron_count_from_density'] - config.explicit_electrons) < 1e-7
        assert abs(result['spin_from_density'] - config.spin) < 1e-7
        np.testing.assert_allclose(coeff.T @ overlap @ coeff, np.eye(coeff.shape[1]), atol=1e-7)
        assert result['ecp_core_electrons'] == config.ecp_core_electrons
        report.append({'element': symbol, 'charge': charge, 'directory': str(directory.relative_to(ROOT)),
                       'energy_hartree': result['total_energy_hartree'], 'converged': result['converged'],
                       'explicit_electrons': result['explicit_electrons'],
                       'electron_trace': result['electron_count_from_density'],
                       'grid_electron_error': result['integration_grid']['electron_count_error'],
                       'scf_cycles': result['scf_cycles'], 'elapsed_seconds': result['elapsed_seconds']})
        (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('REPORT:', output / 'report.json')


if __name__ == '__main__':
    main()
