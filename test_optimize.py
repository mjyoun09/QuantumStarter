"""Geometry optimization: independent-driver agreement, records, refusals, cancellation."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from molecular import OPT_CONVERGENCE, DFTConfig, load_dft_config, parse_coordinates
from scf_control import CalculationCancelled, signature

WATER = 'O 0 0 0; H 0 -.80 .55; H 0 .80 .55'
READY = bool(importlib.util.find_spec('pyscf') and importlib.util.find_spec('geometric'))


class OptimizeInputTests(unittest.TestCase):
    def test_fields(self):
        cfg = DFTConfig(calculation='optimize', opt_convergence='tight', opt_max_steps=20)
        self.assertEqual(DFTConfig.from_dict(json.loads(json.dumps(cfg.to_dict()))), cfg)
        self.assertEqual(signature(cfg), signature(DFTConfig()))
        for kwargs in ({'opt_convergence': 'loose'}, {'opt_max_steps': 0}, {'opt_max_steps': 1.5}):
            with self.subTest(**kwargs), self.assertRaises(ValueError):
                replace(cfg, **kwargs)
        legacy = DFTConfig().to_dict()
        for key in ('opt_convergence', 'opt_max_steps'):
            legacy.pop(key)
        self.assertEqual(DFTConfig.from_dict(legacy), DFTConfig())


@unittest.skipUnless(READY, 'Requires WSL PySCF and geomeTRIC')
class OptimizePhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.cfg = DFTConfig(atoms=parse_coordinates(WATER), xc='pbe', basis='sto-3g', engine='native',
                             calculation='optimize', grid_level=3, cube_points=12, threads=1)

    def tearDown(self):
        self.temp.cleanup()

    def internal(self, xyz):
        r1, r2 = xyz[1] - xyz[0], xyz[2] - xyz[0]
        angle = np.degrees(np.arccos(r1 @ r2 / np.linalg.norm(r1) / np.linalg.norm(r2)))
        return np.array([np.linalg.norm(r1), np.linalg.norm(r2), angle])

    def test_matches_pyscf_geomopt_and_records(self):
        from engines import calculate
        from engines.pyscf_engine import configure_grids
        from pyscf import dft, gto, lib
        from pyscf.geomopt.geometric_solver import optimize
        cfg = replace(self.cfg, opt_convergence='tight')
        result = calculate(cfg, self.root / 'run', lambda _: None)
        opt = result['optimization']
        self.assertEqual(result['status'], 'converged')
        self.assertTrue(opt['converged'])
        final = np.array([a['position'] for a in result['final_geometry_config']['atoms']])
        # Independent driver: PySCF's own SCF, gradient and geomeTRIC interface.
        mol = gto.M(atom=WATER, basis='sto-3g', verbose=0)
        mf = dft.RKS(mol)
        mf.xc, mf.small_rho_cutoff, mf.conv_tol = 'pbe', 0, 1e-10
        configure_grids(mf.grids, cfg)
        c = OPT_CONVERGENCE['tight']
        reference = optimize(mf, maxsteps=100, convergence_energy=c['energy'], convergence_grms=c['grms'],
                             convergence_gmax=c['gmax'], convergence_drms=c['drms'], convergence_dmax=c['dmax'])
        expected = self.internal(reference.atom_coords() * lib.param.BOHR)
        np.testing.assert_allclose(self.internal(final), expected, atol=2e-4)
        # The optimized geometry is a stationary point to the requested force threshold.
        gradient = np.asarray(result['gradient']['gradient_hartree_per_bohr'])
        self.assertLess(np.abs(gradient).max(), c['gmax'])
        self.assertTrue(all(opt['final_gradient_meets_force_criteria'].values()))
        self.assertLess(abs(opt['final_recompute_energy_difference_hartree']), 1e-7)
        energies = [row['energy_hartree'] for row in opt['trajectory']]
        self.assertLess(energies[-1], energies[0])
        run = self.root / 'run'
        self.assertEqual(load_dft_config(run / 'settings.json'), cfg)
        self.assertEqual(result['config'], cfg.to_dict())
        self.assertEqual(load_dft_config(run / 'final_settings.json').calculation, 'gradient')
        frames = (run / 'optimization.xyz').read_text(encoding='utf-8').count('step ')
        self.assertEqual(frames, opt['steps'])
        self.assertTrue((run / 'optimized.xyz').is_file() and (run / 'geometric.log').is_file())
        self.assertEqual(len(list((run / 'steps').iterdir())), opt['steps'])
        second = json.loads((run / 'steps' / 'step_002' / 'results.json').read_text(encoding='utf-8'))
        self.assertEqual(second['initial_guess_details']['resolved'], 'previous_geometry_density')
        self.assertFalse(second['provenance']['restarted'])

    def test_step_limit_scf_failure_cancel_and_missing_library(self):
        from engines import calculate
        limited = calculate(replace(self.cfg, opt_max_steps=1), self.root / 'limited', lambda _: None)
        self.assertEqual(limited['status'], 'not_converged')
        self.assertFalse(limited['converged'])
        self.assertFalse(limited['optimization']['converged'])
        self.assertTrue(limited['warnings'][0].startswith('구조 최적화 미수렴'))
        with self.assertRaisesRegex(RuntimeError, 'SCF가 수렴하지 않아'):
            calculate(replace(self.cfg, max_cycle=1), self.root / 'scf', lambda _: None)
        directory = self.root / 'cancel'

        def cancel_after_first_step(message):
            if message.startswith('OPT_JSON '):
                (directory / 'cancel.request').write_text('', encoding='utf-8')
        with self.assertRaises(CalculationCancelled):
            calculate(self.cfg, directory, cancel_after_first_step)
        self.assertTrue((directory / 'steps' / 'step_001' / 'results.json').is_file())
        with patch.dict(sys.modules, {'geometric': None}):
            with self.assertRaisesRegex(RuntimeError, 'geomeTRIC'):
                calculate(self.cfg, self.root / 'missing', lambda _: None)

    def test_workflow_worker_round_trip(self):
        from dft_workflow import run_calculation
        summary, directory = run_calculation(replace(self.cfg, xc='lda,vwn'), self.root, 'local')
        self.assertEqual(summary['status'], 'converged')
        self.assertEqual(summary['calculation'], 'optimize')
        self.assertTrue((directory / 'optimized.xyz').is_file())


if __name__ == '__main__':
    unittest.main()
