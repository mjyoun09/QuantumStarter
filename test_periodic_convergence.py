"""FFT convergence must not accept a single fortunate interval or failed SCF."""
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from dft_workflow import RunControl
from periodic import load_calculation_config
from periodic_convergence import parse_meshes, run_mesh_convergence, validate_mesh_request


class MeshConvergenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.config = load_calculation_config('examples/silicon_cell.json')
        self.meshes = parse_meshes('21 21 21;25 25 25;29 29 29')

    def run_sweep(self, energies, force_jump=0., failed=None):
        seen = []
        def runner(cfg, root, backend, progress, **kwargs):
            index = len(seen)
            seen.append(cfg)
            directory = Path(root) / str(index)
            directory.mkdir(parents=True)
            forces = [[force_jump * index, 0, 0], [-force_jump * index, 0, 0]]
            result = {'status': 'not_converged' if index == failed else 'converged', 'converged': index != failed,
                      'total_energy_ev': index * 100., 'free_energy_ev': energies[index],
                      'gradient': {'forces_ev_per_angstrom': forces, 'stress_gpa': [[0., 0., 0.]] * 3}}
            return result, directory
        with patch('periodic_convergence.run_calculation', side_effect=runner):
            report, directory = run_mesh_convergence(self.config, self.meshes, self.root, progress=lambda _: None)
        self.assertEqual(json.loads((directory / 'mesh_convergence.json').read_text())['status'], report['status'])
        self.assertTrue((directory / 'mesh_convergence.csv').is_file())
        for cfg in seen:
            self.assertEqual(cfg.atoms, self.config.atoms)
            self.assertEqual(cfg.lattice, self.config.lattice)
            self.assertEqual(cfg.kmesh, self.config.kmesh)
            self.assertEqual(cfg.calculation, 'gradient')
        return report, directory

    def test_two_comparisons_and_free_energy_per_atom(self):
        report, first = self.run_sweep([-2., -2.001, -2.002])
        self.assertTrue(report['converged'])
        self.assertEqual(report['selected_mesh'], [29, 29, 29])
        self.assertAlmostEqual(report['rows'][1]['delta_energy_ev_per_atom'], .0005)
        report, second = self.run_sweep([-1., -2., -2.001])
        self.assertNotEqual(first, second)
        self.assertTrue(first.exists())
        self.assertFalse(report['converged'])
        self.assertIsNone(report['selected_mesh'])
        self.assertTrue(report['rows'][-1]['passed'])

    def test_force_failure_and_scf_failure(self):
        report, _ = self.run_sweep([-2., -2., -2.], force_jump=.1)
        self.assertFalse(report['converged'])
        report, _ = self.run_sweep([-2., -2., -2.], failed=1)
        self.assertEqual(len(report['rows']), 2)
        self.assertFalse(report['converged'])
        self.assertEqual(report['status'], 'not_converged')

    def test_cancellation_and_error_preserve_partial_report(self):
        control = RunControl()
        control.cancel()
        report, directory = run_mesh_convergence(self.config, self.meshes, self.root, control=control)
        self.assertEqual(report['status'], 'cancelled')
        self.assertFalse(report['converged'])
        with patch('periodic_convergence.run_calculation', side_effect=RuntimeError('probe failure')):
            with self.assertRaisesRegex(RuntimeError, 'probe failure'):
                run_mesh_convergence(self.config, self.meshes, self.root, progress=lambda _: None)
        reports = [json.loads(path.read_text()) for path in self.root.glob('*/mesh_convergence.json')]
        self.assertEqual(sorted(r['status'] for r in reports), ['cancelled', 'error'])

    def test_validation(self):
        for meshes in (self.meshes[:2], ((25, 25, 25), (21, 29, 29), (31, 31, 31)),
                       (self.meshes[0], self.meshes[0], self.meshes[2]), ((5, 5), (7, 7), (9, 9)),
                       ((), (7, 7, 7), (9, 9, 9))):
            with self.assertRaises(ValueError):
                validate_mesh_request(self.config, meshes, (.001, .01, .1))
        with self.assertRaises(ValueError):
            validate_mesh_request(self.config, self.meshes, (float('nan'), .01, .1))
        with self.assertRaises(ValueError):
            validate_mesh_request(replace(self.config, calculation='single_point', integral_scheme='gdf', fft_mesh=()),
                                  self.meshes, (.001, .01, .1))


if __name__ == '__main__':
    unittest.main()
