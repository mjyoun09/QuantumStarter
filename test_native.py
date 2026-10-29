"""Native functional derivatives, independently driven SCF and restart contracts."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from molecular import DFTConfig, parse_coordinates
from scf_control import load_checkpoint
from dft_workflow import RunControl, run_calculation, CalculationError, list_dft_runs


class NativeInputTests(unittest.TestCase):
    def test_geometry_and_presets(self):
        from geometry_tools import measure, apply_preset
        atoms = parse_coordinates('H 1 0 0; O 0 0 0; H 0 1 0')
        self.assertAlmostEqual(measure(atoms, [1, 2])[0], 1.)
        self.assertAlmostEqual(measure(atoms, [1, 2], 'bohr')[0], .529177210903)
        self.assertAlmostEqual(measure(atoms, [1, 2, 3])[0], 90.)
        with self.assertRaises(ValueError):
            measure(atoms, [0, 2])
        cfg = DFTConfig()
        preset = apply_preset(cfg, '정밀 격자')
        self.assertEqual(preset.atoms, cfg.atoms)
        self.assertEqual(preset.basis, cfg.basis)
        self.assertEqual(preset.grid_level, 4)

    def test_limits_and_legacy_default(self):
        self.assertEqual(DFTConfig().engine, 'pyscf')
        native = DFTConfig(engine='native')
        self.assertEqual(DFTConfig.from_dict(native.to_dict()), native)
        for kwargs in ({'scf_solver': 'newton'}, {'relativity': 'x2c'}, {'nuclear_model': 'gaussian'},
                       {'native_damping': 1}, {'native_diis_space': 1}, {'native_density_tol': 0},
                       {'atoms': parse_coordinates('Og 0 0 0'), 'basis': 'dyall-v2z'}):
            with self.assertRaises(ValueError):
                replace(native, **kwargs)


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class NativePhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_reference_and_no_external_driver(self):
        from engines import calculate
        from pyscf.dft import rks, uks
        cases = [('He', 'He 0 0 0', 0, 'lda,vwn', 'def2-svp'),
                 ('H2', 'H 0 0 0; H 0 0 .74', 0, 'lda,vwn', 'sto-3g'),
                 ('water', 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587', 0, 'lda,vwn', 'sto-3g'),
                 ('water_pbe', 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587', 0, 'pbe', 'def2-svp'),
                 ('h2_pbe0', 'H 0 0 0; H 0 0 .74', 0, 'pbe0', 'def2-svp'),
                 ('water_b3lyp', 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587', 0, 'b3lyp', 'sto-3g'),
                 ('H', 'H 0 0 0', 1, 'lda,vwn', 'def2-svp'),
                 ('OH', 'O 0 0 0; H 0 0 .97', 1, 'pbe', 'sto-3g'),
                 ('H_pbe0', 'H 0 0 0', 1, 'pbe0', 'def2-svp')]
        for name, atoms, spin, xc, basis in cases:
            with self.subTest(case=name):
                cfg = DFTConfig(atoms=parse_coordinates(atoms), spin=spin, xc=xc, basis=basis,
                                engine='native', grid_level=2, cube_points=12, threads=1, conv_tol=1e-10)
                with patch.object(rks.RKS, 'kernel', side_effect=AssertionError('External SCF called')), \
                     patch.object(uks.UKS, 'kernel', side_effect=AssertionError('External SCF called')):
                    own = calculate(cfg, self.root / name, lambda _: None)
                reference = calculate(replace(cfg, engine='pyscf'), self.root / (name + '_reference'), lambda _: None)
                self.assertTrue(own['converged'], name)
                self.assertTrue(reference['converged'], name)
                if spin == 0:
                    self.assertLess(abs(own['total_energy_hartree'] - reference['total_energy_hartree']), 2e-7)
                else:
                    # Degenerate open-shell orbitals can reach different grid
                    # orientations/SCF states. Compare the same stationary density.
                    import shutil
                    aligned_dir = self.root / (name + '_aligned')
                    aligned_dir.mkdir()
                    shutil.copyfile(self.root / name / 'checkpoint.npz', aligned_dir / 'restart_input.npz')
                    aligned = calculate(replace(cfg, engine='pyscf'), aligned_dir, lambda _: None)
                    self.assertTrue(aligned['converged'])
                    self.assertLess(abs(own['total_energy_hartree'] - aligned['total_energy_hartree']), 2e-7)
                self.assertLess(abs(sum(own['energy_components_hartree'].values()) - own['total_energy_hartree']), 1e-10)
                self.assertLess(own['scf_criteria']['final_residual'], cfg.native_residual_tol)
                self.assertLess(abs(own['scf_history'][-1]['delta_energy_hartree']), cfg.conv_tol)
                self.assertLess(own['scf_history'][-1]['density_change_norm'], cfg.native_density_tol)
                self.assertAlmostEqual(own['scf_history'][-1]['energy_hartree'], own['total_energy_hartree'], places=10)
                self.assertAlmostEqual(own['electron_count_from_density'], cfg.electrons, places=8)
                self.assertAlmostEqual(own['spin_from_density'], cfg.spin, places=8)
                self.assertEqual(own['provenance']['scf_driver'], 'QuantumStarter')
                data = np.load(self.root / name / 'electronic_structure.npz')
                c, s, occ, dm = (data[k] for k in ('mo_coeff', 'overlap', 'mo_occupation', 'density_matrix'))
                np.testing.assert_allclose(c.swapaxes(-1, -2) @ s @ c, np.broadcast_to(np.eye(c.shape[-1]), dm.shape), atol=1e-8)
                np.testing.assert_allclose((c * occ[..., None, :]) @ c.swapaxes(-1, -2), dm, atol=1e-10)

    def test_fock_is_energy_derivative_rks_and_uks(self):
        from pyscf import gto, lib
        from engines.native_scf import NativeKS
        from engines.pyscf_engine import configure_grids
        lib.num_threads(1)
        for method in ('rks', 'uks'):
            for xc in ('lda,vwn', 'pbe', 'pbe0', 'b3lyp'):
                with self.subTest(method=method, xc=xc):
                    cfg = DFTConfig(atoms=parse_coordinates('H 0 0 0; H 0 0 .74'), basis='sto-3g',
                                    engine='native', method=method, xc=xc, grid_level=2)
                    mol = gto.M(atom='H 0 0 0; H 0 0 .74', basis=cfg.basis, verbose=0)
                    solver = NativeKS(mol, cfg, self.root)
                    configure_grids(solver.grids, cfg)
                    solver.grids.build()
                    dm = np.linalg.inv(solver.overlap)
                    perturb = np.array([[.2, .1], [.1, -.15]])
                    if method == 'uks':
                        dm = np.stack((dm * .6, dm * .4))
                        perturb = np.stack((perturb, perturb * -.3))
                    _, fock, _ = solver.evaluate(dm)
                    step = 1e-5
                    derivative = (solver.evaluate(dm + step * perturb)[0] - solver.evaluate(dm - step * perturb)[0]) / (2 * step)
                    self.assertAlmostEqual(derivative, float(np.sum(fock * perturb)), delta=2e-7)

    def test_guesses_grid_basis_and_nonconvergence(self):
        from engines import calculate
        cfg = DFTConfig(engine='native', basis='sto-3g', xc='lda,vwn', grid_level=2, threads=1, cube_points=12)
        first = calculate(cfg, self.root / 'first', lambda _: None)
        second = calculate(replace(cfg, initial_guess='perturbed'), self.root / 'perturbed', lambda _: None)
        self.assertTrue(second['converged'])
        self.assertAlmostEqual(first['total_energy_hartree'], second['total_energy_hartree'], places=8)
        fine = calculate(replace(cfg, grid_level=3), self.root / 'fine', lambda _: None)
        self.assertTrue(fine['converged'])
        self.assertLess(abs(fine['total_energy_hartree'] - first['total_energy_hartree']), 1e-4)
        failed = calculate(replace(cfg, max_cycle=1), self.root / 'incomplete', lambda _: None)
        self.assertEqual(failed['status'], 'not_converged')
        self.assertTrue((self.root / 'incomplete' / 'checkpoint.npz').exists())

    def test_cancel_restart_and_mismatch(self):
        cfg = DFTConfig(engine='native', basis='sto-3g', xc='lda,vwn', grid_level=1, cube_points=12, threads=1)
        control = RunControl()
        def progress(text):
            if text.startswith('SCF_JSON '):
                control.cancel()
        result, directory = run_calculation(cfg, self.root, 'local', progress, control=control)
        self.assertEqual(result['status'], 'cancelled')
        self.assertFalse(result['converged'])
        self.assertTrue((directory / 'checkpoint.npz').exists())
        self.assertEqual(list_dft_runs(self.root)[0]['summary']['status'], 'cancelled')
        restarted, final = run_calculation(cfg, self.root, 'local', restart_from=directory)
        self.assertTrue(restarted['converged'])
        self.assertTrue(restarted['provenance']['restarted'])
        self.assertTrue((final / 'restart_input.npz').exists())
        with self.assertRaises(CalculationError):
            run_calculation(replace(cfg, xc='pbe'), self.root, 'local', restart_from=directory)

    def test_near_degenerate_open_shell_converges(self):
        # Regression (2026-10-29): a coefficient cap of 30 reset the Pulay history on
        # ~289 of 300 cycles for this OH case, so the SCF stalled at dD ~2e-5.
        from engines import calculate
        from engines.native_scf import PulayDIIS
        diis = PulayDIIS(8)
        base = np.eye(2)
        for k in range(4):
            out = diis.update(base * (1 + k), np.array([[1e-6, 1e-6 * (1 + 1e-9 * k)], [0, 0]]))
        self.assertTrue(np.isfinite(out).all())
        self.assertEqual(len(diis.focks), 4)
        cfg = DFTConfig(atoms=parse_coordinates('O 0 0 0; H 0.1 0 .97'), spin=1, xc='pbe', basis='sto-3g',
                        engine='native', grid_level=1, conv_tol=1e-12, native_density_tol=1e-10,
                        native_residual_tol=1e-9, max_cycle=100, cube_points=12, threads=1)
        result = calculate(cfg, self.root / 'oh', lambda _: None)
        self.assertTrue(result['converged'])
        self.assertLess(result['scf_cycles'], 60)
        # Converged means stationary, not lowest: this OH state is internally unstable
        # (PySCF stability analysis finds a lower UKS solution ~8e-6 Ha below).

    def test_atomic_density_guess_and_late_damping(self):
        # 0.8.6: auto guess = superposition of atomic densities; mixing only while the
        # residual is large. Water dropped from 14 to 9 cycles; the Fe atom (2S=4),
        # which never converged with the core guess, now converges.
        from engines import calculate
        water = DFTConfig(atoms=parse_coordinates('O 0 0 0; H 0 -.757 .587; H 0 .757 .587'), xc='pbe',
                          basis='def2-svp', engine='native', grid_level=3, cube_points=12, threads=1)
        result = calculate(water, self.root / 'water', lambda _: None)
        self.assertEqual(result['initial_guess_details']['resolved'], 'minao')
        self.assertLessEqual(result['scf_cycles'], 11)
        core = calculate(replace(water, initial_guess='core'), self.root / 'core', lambda _: None)
        self.assertAlmostEqual(core['total_energy_hartree'], result['total_energy_hartree'], delta=1e-9)
        iron = calculate(DFTConfig(atoms=parse_coordinates('Fe 0 0 0'), spin=4, xc='pbe', basis='sto-3g',
                                   engine='native', grid_level=3, cube_points=12, threads=1, max_cycle=250),
                         self.root / 'fe', lambda _: None)
        self.assertTrue(iron['converged'])

    def test_external_checkpoint_and_restart(self):
        cfg = DFTConfig(basis='sto-3g', xc='lda,vwn', grid_level=1, cube_points=12, threads=1, max_cycle=1)
        for solver in ('diis', 'newton'):
            first, directory = run_calculation(replace(cfg, scf_solver=solver), self.root, 'local')
            self.assertTrue((directory / 'checkpoint.npz').exists())
            final, _ = run_calculation(replace(cfg, scf_solver=solver, max_cycle=100), self.root, 'local', restart_from=directory)
            self.assertTrue(final['converged'])
            self.assertTrue(final['provenance']['restarted'])


if __name__ == '__main__':
    unittest.main()
