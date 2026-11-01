"""Cell-coordinate/force/stress consistency, termination and failure contracts."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms
from periodic import PeriodicConfig, load_calculation_config
from scf_control import CalculationCancelled


class CellInputTests(unittest.TestCase):
    def test_validation_and_roundtrip(self):
        config = load_calculation_config('examples/silicon_cell.json')
        self.assertEqual(PeriodicConfig.from_dict(config.to_dict()), config)
        for changes in ({'cell_max_steps': 0}, {'cell_relaxation': 'unknown'},
                        {'cell_force_tol': float('nan')}, {'cell_energy_tol_ev': -1},
                        {'cell_stress_tol_gpa': 0}, {'integral_scheme': 'gdf', 'fft_mesh': ()}):
            with self.subTest(**changes), self.assertRaises(ValueError):
                replace(config, **changes)


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF import')
class CellOptimizerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.target = np.eye(3) * 4
        self.fractional = np.array([[.25, .25, .25]])
        self.seen = []

    def config(self, mode='full', **changes):
        deform = np.array([[1.05, .02, 0], [.02, .97, 0], [0, 0, 1.02]]) if mode == 'full' else np.eye(3) * 1.05
        cell = self.target @ deform
        cfg = PeriodicConfig.from_ase(Atoms('Si', positions=self.fractional @ cell + [.04, -.02, 0], cell=cell, pbc=True),
                    integral_scheme='fftdf', fft_mesh=(17, 19, 21), kmesh=(2, 1, 1),
                    calculation='cell_optimize', cell_relaxation=mode, cell_max_steps=80,
                    cell_force_tol=1e-4, cell_stress_tol_gpa=1e-3, cell_energy_tol_ev=1e-8, smearing_ev=.1)
        return replace(cfg, **changes)

    def elastic(self, cfg, directory, progress, **kwargs):
        """Independent harmonic cell + internal-coordinate potential, analytic derivatives."""
        from engines.periodic_optimizer import EV_PER_ANGSTROM3_GPA
        self.seen.append(cfg)
        cell, positions = np.asarray(cfg.lattice), np.asarray([a.position for a in cfg.atoms])
        displacement = positions - self.fractional @ cell
        delta = cell - self.target
        energy = float(2.5 * np.sum(delta ** 2) + 1.5 * np.sum(displacement ** 2))
        forces = -3 * displacement
        virial = 5 * cell.T @ delta + 3 * displacement.T @ displacement
        stress = (virial + virial.T) / (2 * np.linalg.det(cell)) * EV_PER_ANGSTROM3_GPA
        Path(directory).mkdir(parents=True, exist_ok=True)
        return {'converged': True, 'status': 'converged', 'config': cfg.to_dict(), 'warnings': [],
                'free_energy_ev': energy, 'total_energy_ev': energy + np.linalg.det(cell),
                'gradient': {'forces_ev_per_angstrom': forces.tolist(), 'stress_gpa': stress.tolist(),
                             'max_atomic_force_ev_per_angstrom': float(np.linalg.norm(forces, axis=1).max())}}

    def run_elastic(self, config, label='run'):
        from engines.periodic_optimizer import optimize_periodic
        with patch('engines.periodic_optimizer.calculate_periodic', side_effect=self.elastic):
            return optimize_periodic(config, self.root / label, lambda _: None)

    def test_full_cell_shear_and_atomic_relaxation(self):
        result = self.run_elastic(self.config())
        self.assertTrue(result['converged'])
        self.assertLess(result['free_energy_ev'], 1e-8)
        np.testing.assert_allclose(result['final_geometry_config']['lattice'], self.target, atol=2e-5)
        self.assertTrue(all(c.fft_mesh == (17, 19, 21) and c.kmesh == (2, 1, 1) for c in self.seen))
        self.assertTrue(all(c.calculation == 'gradient' for c in self.seen))
        self.assertEqual(result['config']['calculation'], 'cell_optimize')
        for filename in ('settings.json', 'final_settings.json', 'results.json', 'cell_history.json',
                         'cell_optimization.csv', 'cell_optimization.extxyz', 'optimized.extxyz'):
            self.assertTrue((self.root / 'run' / filename).is_file(), filename)
        block = result['cell_optimization']
        self.assertIn('free energy F', block['objective'])
        self.assertLess(block['final_objective_ev'], block['initial_objective_ev'])

    def test_isotropic_keeps_shape(self):
        cfg = self.config('isotropic')
        result = self.run_elastic(cfg)
        self.assertTrue(result['converged'])
        final = np.asarray(result['final_geometry_config']['lattice'])
        np.testing.assert_allclose(final / np.linalg.det(final) ** (1 / 3), np.eye(3), atol=1e-12)

    def test_automatic_mesh_is_resolved_once(self):
        cfg = self.config('isotropic', fft_mesh=())
        with patch('engines.periodic_optimizer.build_cell', return_value=(SimpleNamespace(mesh=[17, 19, 21]), {})) as build:
            result = self.run_elastic(cfg)
        build.assert_called_once_with(cfg)
        self.assertTrue(result['converged'])
        self.assertEqual(result['config']['fft_mesh'], ())
        self.assertTrue(all(c.fft_mesh == (17, 19, 21) for c in self.seen))

    def test_step_limit_and_final_recheck(self):
        limited = self.run_elastic(self.config(cell_max_steps=1), 'limited')
        self.assertFalse(limited['converged'])
        self.assertEqual(limited['cell_optimization']['steps'], 1)
        from engines.periodic_optimizer import optimize_periodic
        def changed_final(cfg, directory, progress, **kwargs):
            result = self.elastic(cfg, directory, progress, **kwargs)
            if Path(directory).name == 'final':
                result['free_energy_ev'] += .5
            return result
        with patch('engines.periodic_optimizer.calculate_periodic', side_effect=changed_final):
            result = optimize_periodic(self.config(), self.root / 'changed', lambda _: None)
        self.assertFalse(result['converged'])
        self.assertFalse(result['cell_optimization']['final_criteria_met']['energy_change'])

    def test_nonconverged_scf_and_cancellation_are_not_success(self):
        from engines.periodic_optimizer import optimize_periodic
        with patch('engines.periodic_optimizer.calculate_periodic', return_value={'converged': False, 'gradient': None}):
            with self.assertRaisesRegex(RuntimeError, 'SCF 미수렴'):
                optimize_periodic(self.config(), self.root / 'failed', lambda _: None)
        self.assertFalse((self.root / 'failed' / 'results.json').exists())
        target = self.root / 'cancelled'
        target.mkdir()
        (target / 'cancel.request').touch()
        with self.assertRaises(CalculationCancelled):
            optimize_periodic(self.config(), target, lambda _: None)

    def test_stress_residual_mode(self):
        from engines.periodic_optimizer import stress_residual
        stress = np.diag([1., -1., 0.])
        self.assertEqual(stress_residual(stress, 'isotropic'), 0)
        self.assertEqual(stress_residual(stress, 'full'), 1)

    def test_final_cancel_does_not_publish_single_point_success(self):
        from engines.periodic_optimizer import optimize_periodic
        cfg = self.config()
        directory = self.root / 'final_cancel'
        def cancelled_final(config, path, progress, **kwargs):
            result = self.elastic(config, path, progress, **kwargs)
            if Path(path).name == 'final':
                (directory / 'cancel.request').touch()
            return result
        with patch('engines.periodic_optimizer.calculate_periodic', side_effect=cancelled_final):
            with self.assertRaises(CalculationCancelled):
                optimize_periodic(cfg, directory, lambda _: None)
        self.assertFalse((directory / 'results.json').exists())
        self.assertEqual(load_calculation_config(directory / 'settings.json'), cfg)
        self.assertTrue((directory / 'cell_history.json').is_file())


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class CellWorkerTests(unittest.TestCase):
    def test_cancel_between_steps_and_scf_failure(self):
        from dft_workflow import CalculationError, RunControl, run_calculation
        cfg = replace(load_calculation_config('examples/silicon_cell.json'), kmesh=(1, 1, 1), fft_mesh=(17, 17, 17))
        with tempfile.TemporaryDirectory() as temp:
            control = RunControl()
            def progress(message):
                if message.startswith('셀 최적화 0:'):
                    control.cancel()
            summary, directory = run_calculation(cfg, temp, progress=progress, control=control)
            self.assertEqual(summary['status'], 'cancelled')
            self.assertEqual(summary['config']['calculation'], 'cell_optimize')
            self.assertTrue((directory / 'cell_history.json').is_file())
            self.assertFalse((directory / 'final_settings.json').exists())
            with self.assertRaises(CalculationError) as caught:
                run_calculation(replace(cfg, max_cycle=1), temp)
            self.assertTrue((caught.exception.directory / 'error.json').is_file())
            failed = json.loads((caught.exception.directory / 'steps/step_0000/results.json').read_text())
            self.assertFalse(failed['converged'])
            self.assertIsNone(failed['gradient'])


class CellCliTests(unittest.TestCase):
    def test_cell_options_reach_worker(self):
        import main
        args = ['main.py', '--config', 'examples/silicon_cell.json', '--calculation', 'cell_optimize',
                '--cell-relaxation', 'full', '--cell-max-steps', '2', '--cell-force-tol', '.03',
                '--cell-stress-tol-gpa', '.2', '--cell-energy-tol-ev', '1e-6', '--no-show']
        with patch('sys.argv', args), patch('dft_workflow.run_calculation',
                                           return_value=({'status': 'cancelled'}, Path('unused'))) as runner:
            self.assertEqual(main.main(), 3)
        cfg = runner.call_args.args[0]
        self.assertEqual(cfg.cell_max_steps, 2)
        self.assertEqual(cfg.cell_relaxation, 'full')
        self.assertEqual(cfg.cell_force_tol, .03)
        self.assertEqual(cfg.cell_stress_tol_gpa, .2)
        self.assertEqual(cfg.cell_energy_tol_ev, 1e-6)
        for args in (['main.py', '--cell-max-steps', '2'],
                     ['main.py', '--hydrogen-check', '--fft-mesh', '25', '25', '25']):
            with self.subTest(args=args), patch('sys.argv', args), patch('sys.stderr'), self.assertRaises(SystemExit) as caught:
                main.main()
            self.assertEqual(caught.exception.code, 2)


if __name__ == '__main__':
    unittest.main()
