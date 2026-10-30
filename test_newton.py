"""Native trust-region Newton SCF: same minimum as DIIS, escape from an unstable (saddle) solution."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from molecular import DFTConfig, parse_coordinates


class TrustRegionTests(unittest.TestCase):
    def test_zero_gradient_and_stiff_positive_quadratic(self):
        from engines.native_newton import steihaug
        diagonal = np.array([1e-4, 1000.])
        gradient = np.array([1e-5, .1])
        step = steihaug(lambda x: diagonal * x, gradient, 1., 1e-12, diagonal=diagonal)
        np.testing.assert_allclose(step, -gradient / diagonal, atol=1e-12)
        np.testing.assert_array_equal(steihaug(lambda x: x, np.zeros(0), 1., 0.), np.zeros(0))

    def test_negative_curvature_reaches_boundary_downhill(self):
        from engines.native_newton import steihaug
        diagonal, gradient = np.array([-2., 3.]), np.array([.1, 0.])
        step = steihaug(lambda x: diagonal * x, gradient, .2, 1e-12)
        self.assertAlmostEqual(np.linalg.norm(step), .2)
        self.assertLess(gradient @ step + .5 * step @ (diagonal * step), 0.)


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class NativeNewtonTests(unittest.TestCase):
    def setUp(self):
        from pyscf import lib
        old_threads = lib.num_threads()
        lib.num_threads(1)
        self.addCleanup(lib.num_threads, old_threads)
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_water_matches_diis(self):
        from engines import calculate
        cfg = DFTConfig(atoms=parse_coordinates('O 0 0 0; H 0 -.757 .587; H 0 .757 .587'), xc='pbe', basis='sto-3g',
                        engine='native', scf_solver='newton', grid_level=2, cube_points=12, threads=1)
        newton = calculate(cfg, self.root / 'newton', lambda _: None)
        diis = calculate(replace(cfg, scf_solver='diis'), self.root / 'diis', lambda _: None)
        self.assertTrue(newton['converged'])
        self.assertAlmostEqual(newton['total_energy_hartree'], diis['total_energy_hartree'], delta=1e-8)
        self.assertLess(newton['scf_criteria']['final_residual'], cfg.native_residual_tol)
        self.assertLessEqual(newton['native_newton']['aufbau_violation_hartree'], 0)
        self.assertAlmostEqual(newton['electron_count_from_density'], 10, places=8)

    def test_escapes_unstable_solution(self):
        # DIIS converges OH (coarse grid) to a state that the stability analysis shows to be a
        # saddle point. Started there, Newton must leave it and stop only at a local minimum.
        from pyscf import gto
        from engines.native_scf import NativeKS
        from engines.native_stability import OrbitalRotations, lowest_mode
        from engines.pyscf_engine import configure_grids
        cfg = DFTConfig(atoms=parse_coordinates('O 0 0 0; H 0.1 0 .97'), spin=1, xc='pbe', basis='sto-3g',
                        engine='native', grid_level=1, conv_tol=1e-12, native_density_tol=1e-10,
                        native_residual_tol=1e-9, max_cycle=300)
        mol = gto.M(atom='O 0 0 0; H 0.1 0 .97', basis='sto-3g', spin=1, verbose=0)
        diis = NativeKS(mol, cfg, self.root)
        configure_grids(diis.grids, cfg)
        saddle = diis.kernel()
        self.assertLess(lowest_mode(OrbitalRotations(diis))[0], -1e-4)
        newton = NativeKS(mol, replace(cfg, scf_solver='newton', max_cycle=100), self.root)
        configure_grids(newton.grids, cfg)
        energy = newton.kernel(dm0=diis.make_rdm1())
        self.assertTrue(newton.converged)
        self.assertLess(energy, saddle - 1e-7)
        # Negative curvature can be followed by the trust-region step itself;
        # an explicit post-convergence escape is not required on every trajectory.
        self.assertGreater(newton.newton_escapes[-1]['lowest_eigenvalue_hartree'], -1e-4)
        self.assertGreater(lowest_mode(OrbitalRotations(newton))[0], -1e-4)
        self.assertAlmostEqual(np.trace(newton.make_rdm1().sum(axis=0) @ newton.overlap), 9, places=8)

    def test_transition_metals_converge_to_stable_integer_occupations(self):
        from pyscf import dft, gto
        from engines.native_scf import NativeKS
        from engines.pyscf_engine import configure_grids
        for element, spin in [('Cr', 6), ('Co', 3)]:
            with self.subTest(element=element):
                cfg = DFTConfig(atoms=parse_coordinates(f'{element} 0 0 0'), spin=spin, xc='pbe',
                                basis='sto-3g', engine='native', scf_solver='newton', grid_level=3,
                                max_cycle=100, threads=1)
                mol = gto.M(atom=f'{element} 0 0 0', spin=spin, basis=cfg.basis, verbose=0)
                mf = NativeKS(mol, cfg, self.root)
                configure_grids(mf.grids, cfg)
                energy = mf.kernel()
                self.assertTrue(mf.converged)
                self.assertLess(mf.final_residual, cfg.native_residual_tol)
                self.assertGreater(mf.newton_escapes[-1]['lowest_eigenvalue_hartree'], -1e-4)
                self.assertTrue(np.all((mf.mo_occ == 0) | (mf.mo_occ == 1)))
                np.testing.assert_allclose(np.einsum('sij,ji->s', mf.make_rdm1(), mf.overlap),
                                           mol.nelec, atol=1e-9)
                reference = dft.UKS(mol)
                reference.xc, reference.grids, reference.small_rho_cutoff = cfg.xc, mf.grids, 0
                reference = reference.newton()
                reference.conv_tol, reference.conv_tol_grad = 1e-11, 1e-8
                ref_energy = reference.kernel(dm0=mf.make_rdm1())
                self.assertTrue(reference.converged)
                self.assertAlmostEqual(energy, ref_energy, delta=1e-9)

    def test_filled_basis_and_forced_nonconvergence(self):
        from engines import calculate
        cfg = DFTConfig(atoms=parse_coordinates('He 0 0 0'), xc='pbe', basis='sto-3g', engine='native',
                        scf_solver='newton', cube_points=12, threads=1)
        result = calculate(cfg, self.root / 'filled', lambda _: None)
        self.assertTrue(result['converged'])
        self.assertEqual(result['native_newton']['saddle_escapes'][-1]['lowest_eigenvalue_hartree'], 0)
        self.assertIn('analytic', result['native_newton']['hessian_vector'])
        failed = calculate(replace(cfg, max_cycle=1), self.root / 'one', lambda _: None)
        self.assertFalse(failed['converged'])
        one_rotation = calculate(replace(cfg, atoms=parse_coordinates('H 0 0 0; H 0 0 .74')),
                                 self.root / 'one_rotation', lambda _: None)
        self.assertTrue(one_rotation['converged'])
        self.assertGreater(one_rotation['native_newton']['saddle_escapes'][-1]['lowest_eigenvalue_hartree'], 0)


if __name__ == '__main__':
    unittest.main()
