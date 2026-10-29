"""Nuclear gradients: frozen-grid finite differences, engine agreement, D4, ECP, refusals."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from molecular import DFTConfig, parse_coordinates
from scf_control import signature

WATER = 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587'
# Bent doublet radical with a non-degenerate singly occupied orbital. (OH, whose pi pair is
# nearly degenerate, oscillates in the native DIIS on coarse grids; that is an SCF issue.)
NH2 = 'N 0 0 0; H 0 -.80 .62; H 0.05 .80 .60'
HAS_PYSCF = bool(importlib.util.find_spec('pyscf'))


class GradientInputTests(unittest.TestCase):
    def test_calculation_choices(self):
        cfg = DFTConfig(calculation='gradient')
        self.assertEqual(DFTConfig.from_dict(cfg.to_dict()), cfg)
        self.assertEqual(signature(cfg), signature(DFTConfig()))
        for value in ('scan', 'Gradient', ''):
            with self.subTest(value=value), self.assertRaises(ValueError):
                DFTConfig(calculation=value)


@unittest.skipUnless(HAS_PYSCF, 'Requires WSL PySCF')
class GradientPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def native(self, mol, cfg, grids=None, dm0=None):
        from engines.native_scf import NativeKS
        from engines.pyscf_engine import configure_grids
        mf = NativeKS(mol, cfg, self.root)
        if grids is None:
            configure_grids(mf.grids, cfg)
        else:
            mf.grids = grids  # Bound to the reference molecule: rebuilding gives the same points.
        energy = mf.kernel(dm0=dm0)
        self.assertTrue(mf.converged)
        return mf, energy

    def test_native_is_exact_derivative_on_frozen_grid(self):
        from pyscf import gto, lib
        from engines.native_gradient import native_gradient
        lib.num_threads(1)
        step = 1e-4
        cases = [(WATER, 0, 'lda,vwn', [(0, 2), (1, 1), (1, 2), (2, 0)]),
                 (WATER, 0, 'pbe0', [(0, 2), (1, 1), (2, 2)]),
                 (NH2, 1, 'pbe', [(0, 0), (0, 2), (2, 0), (2, 1)]),
                 (NH2, 1, 'b3lyp', [(0, 1), (1, 2)])]
        for atoms, spin, xc, coordinates in cases:
            cfg = DFTConfig(atoms=parse_coordinates(atoms), spin=spin, xc=xc, basis='sto-3g', engine='native',
                            grid_level=1, conv_tol=1e-12, native_density_tol=1e-10, native_residual_tol=1e-9,
                            max_cycle=300)
            mol = gto.M(atom=atoms, basis='sto-3g', spin=spin, verbose=0)
            reference, _ = self.native(mol, cfg)
            analytic = sum(native_gradient(reference).values())
            for atom, axis in coordinates:
                with self.subTest(xc=xc, spin=spin, atom=atom, axis=axis):
                    energies = []
                    for sign in (1, -1):
                        coords = mol.atom_coords().copy()
                        coords[atom, axis] += sign * step
                        moved = mol.set_geom_(coords, unit='Bohr', inplace=False)
                        energies.append(self.native(moved, cfg, reference.grids, reference.make_rdm1())[1])
                    numeric = (energies[0] - energies[1]) / (2 * step)
                    self.assertAlmostEqual(analytic[atom, axis], numeric, delta=1e-6)

    def test_engines_agree_and_result_contract(self):
        from engines import calculate
        cfg = DFTConfig(atoms=parse_coordinates(WATER), basis='def2-svp', xc='pbe', calculation='gradient',
                        grid_level=2, cube_points=12, threads=1, conv_tol=1e-11,
                        native_density_tol=1e-9, native_residual_tol=1e-8, max_cycle=200)
        external = calculate(cfg, self.root / 'pyscf', lambda _: None)
        own = calculate(replace(cfg, engine='native'), self.root / 'native', lambda _: None)
        g_ext = np.asarray(external['gradient']['gradient_hartree_per_bohr'])
        g_own = np.asarray(own['gradient']['gradient_hartree_per_bohr'])
        self.assertEqual((external['gradient']['driver'], own['gradient']['driver']), ('PySCF', 'QuantumStarter'))
        self.assertLess(np.abs(g_ext - g_own).max(), 1e-6)
        for result, gradient in ((external, g_ext), (own, g_own)):
            block = result['gradient']
            components = sum(np.asarray(v) for v in block['components_hartree_per_bohr'].values())
            np.testing.assert_allclose(components, gradient, atol=1e-14)
            forces = np.asarray(block['forces_ev_per_angstrom'])
            np.testing.assert_allclose(forces, -gradient * 27.211386245988 / 0.529177210903, rtol=1e-9)
            self.assertAlmostEqual(block['max_atomic_force_ev_per_angstrom'],
                                   np.linalg.norm(forces, axis=1).max(), places=12)
            self.assertLess(np.abs(block['net_force_hartree_per_bohr']).max(), 1e-3)
            self.assertFalse(block['grid_response'])
            # Mirror symmetry of water (x = 0 plane, H1 <-> H2 under y -> -y).
            self.assertLess(np.abs(gradient[:, 0]).max(), 1e-8)
            self.assertAlmostEqual(gradient[1, 1], -gradient[2, 1], delta=1e-8)
        self.assertEqual(set(own['gradient']['components_hartree_per_bohr']),
                         {'one_electron', 'coulomb', 'exact_exchange', 'xc', 'overlap_pulay', 'nuclear_repulsion'})
        single = calculate(replace(cfg, calculation='single_point'), self.root / 'single', lambda _: None)
        self.assertIsNone(single['gradient'])
        self.assertAlmostEqual(single['total_energy_hartree'], external['total_energy_hartree'], delta=1e-10)

    def test_moving_grid_difference_is_small(self):
        # The physical energy moves its grid with the atoms; the analytic result omits that
        # (grid-response) term. Its size here bounds the approximation, it is not zero.
        from pyscf import gto
        from engines.native_gradient import native_gradient
        cfg = DFTConfig(atoms=parse_coordinates(WATER), basis='sto-3g', xc='pbe', engine='native', grid_level=3,
                        conv_tol=1e-12, native_density_tol=1e-10, native_residual_tol=1e-9, max_cycle=300)
        mol = gto.M(atom=WATER, basis='sto-3g', verbose=0)
        reference, _ = self.native(mol, cfg)
        analytic = sum(native_gradient(reference).values())
        energies = []
        for sign in (1, -1):
            coords = mol.atom_coords().copy()
            coords[0, 2] += sign * 1e-3
            energies.append(self.native(mol.set_geom_(coords, unit='Bohr', inplace=False), cfg,
                                        dm0=reference.make_rdm1())[1])
        self.assertLess(abs(analytic[0, 2] - (energies[0] - energies[1]) / 2e-3), 1e-4)

    @unittest.skipUnless(importlib.util.find_spec('dftd4'), 'Requires dftd4')
    def test_d4_gradient_and_total(self):
        from engines import calculate
        from engines.dispersion import d4_energy
        numbers = [8, 1, 1, 18]
        coords = np.array([[0., 0, 0], [0, -1.43, 1.11], [0, 1.43, 1.11], [0.4, 0.3, 6.5]])
        analytic = d4_energy(numbers, coords, 0, 'pbe', gradient=True)['gradient_hartree_per_bohr']
        step = 1e-5
        for atom in range(4):
            for axis in range(3):
                plus, minus = coords.copy(), coords.copy()
                plus[atom, axis] += step
                minus[atom, axis] -= step
                numeric = (d4_energy(numbers, plus, 0, 'pbe')['energy_hartree']
                           - d4_energy(numbers, minus, 0, 'pbe')['energy_hartree']) / (2 * step)
                self.assertAlmostEqual(analytic[atom, axis], numeric, delta=1e-9)
        cfg = DFTConfig(atoms=parse_coordinates(WATER), basis='sto-3g', xc='pbe', calculation='gradient',
                        dispersion='d4', grid_level=1, cube_points=12, threads=1)
        for engine in ('pyscf', 'native'):
            result = calculate(replace(cfg, engine=engine), self.root / engine, lambda _: None)
            block = result['gradient']
            self.assertTrue(block['includes_dispersion'])
            parts = block['components_hartree_per_bohr']
            bohr = np.array([a.position for a in cfg.atoms]) / 0.52917721092
            expected = d4_energy([8, 1, 1], bohr, 0, 'pbe', gradient=True)['gradient_hartree_per_bohr']
            np.testing.assert_allclose(parts['dispersion'], expected, atol=1e-12)

    def test_ecp_path_and_unconverged_refusal(self):
        from engines import calculate
        auh = DFTConfig(atoms=parse_coordinates('Au 0 0 0; H 0 0 1.52'), basis='def2-svp', xc='pbe',
                        calculation='gradient', grid_level=2, cube_points=12, threads=1)
        result = calculate(auh, self.root / 'auh', lambda _: None)
        gradient = np.asarray(result['gradient']['gradient_hartree_per_bohr'])
        self.assertTrue(result['converged'])
        self.assertEqual(result['ecp']['core_electrons_by_atom'], [60, 0])
        self.assertLess(np.abs(gradient[:, :2]).max(), 1e-8)
        self.assertLess(abs(gradient[0, 2] + gradient[1, 2]), 1e-3)
        failed = calculate(replace(auh, atoms=parse_coordinates(WATER), basis='sto-3g', max_cycle=1),
                           self.root / 'failed', lambda _: None)
        self.assertEqual(failed['status'], 'not_converged')
        self.assertIsNone(failed['gradient'])
        self.assertTrue(any('힘을 계산하지 않았습니다' in w for w in failed['warnings']))


if __name__ == '__main__':
    unittest.main()
