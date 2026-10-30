"""Bulk FFTDF forces and stress through the product API: finite differences, units, smearing, mixing."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from periodic import load_calculation_config

HARTREE_EV, BOHR_A = 27.211386245988, 0.529177210903


def silicon(**settings):
    cfg = load_calculation_config('examples/silicon_bulk.json')
    cfg = replace(cfg, atoms=(cfg.atoms[0], replace(cfg.atoms[1], position=(1.4375, 1.3575, 1.3575))))
    base = dict(basis='gth-szv-molopt-sr', kmesh=(1, 1, 1), precision=1e-8, conv_tol=1e-11, smearing_ev=0.0,
                integral_scheme='fftdf', fft_mesh=(25, 25, 25), threads=1)
    base.update(settings)
    return replace(cfg, **base)


class PeriodicGradientInputTests(unittest.TestCase):
    def test_rules(self):
        cfg = silicon(calculation='gradient')
        self.assertEqual((cfg.integral_scheme, cfg.fft_mesh), ('fftdf', (25, 25, 25)))
        for kwargs, message in (({'integral_scheme': 'gdf', 'fft_mesh': ()}, 'fftdf'),
                                ({'integral_scheme': 'gdf', 'calculation': 'single_point'}, 'fft_mesh'),
                                ({'fft_mesh': (2, 25, 25)}, 'FFT mesh'), ({'calculation': 'optimize'}, 'calculation')):
            with self.subTest(**kwargs), self.assertRaisesRegex(ValueError, message):
                replace(cfg, **kwargs)


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class PeriodicGradientPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def run_cfg(self, cfg, name):
        from engines import calculate
        result = calculate(cfg, self.root / name, lambda _: None)
        self.assertTrue(result['converged'], name)
        return result

    def displaced(self, cfg, step):
        atoms = list(cfg.atoms)
        x, y, z = atoms[1].position
        atoms[1] = replace(atoms[1], position=(x + step * BOHR_A, y, z))
        return replace(cfg, atoms=tuple(atoms), calculation='single_point')

    def strained(self, cfg, step):
        strain = np.eye(3)
        strain[0, 0] += step
        lattice = np.asarray(cfg.lattice) @ strain
        atoms = tuple(replace(a, position=tuple(np.asarray(a.position) @ strain)) for a in cfg.atoms)
        return replace(cfg, lattice=tuple(map(tuple, lattice)), atoms=atoms, calculation='single_point')

    def test_forces_and_stress_match_finite_differences(self):
        cfg = silicon(calculation='gradient')
        result = self.run_cfg(cfg, 'ref')
        block = result['gradient']
        self.assertEqual(block['derivative_of'], 'total energy E')
        self.assertEqual(result['integral_scheme'], 'FFTDF')
        self.assertEqual(result['fft_mesh'], [25, 25, 25])
        h = 1e-3  # bohr
        numeric = (self.run_cfg(self.displaced(cfg, h), 'p')['total_energy_hartree']
                   - self.run_cfg(self.displaced(cfg, -h), 'm')['total_energy_hartree']) / (2 * h)
        self.assertAlmostEqual(block['gradient_hartree_per_bohr'][1][0], numeric, delta=1e-7)
        forces = np.asarray(block['forces_ev_per_angstrom'])
        np.testing.assert_allclose(forces, -np.asarray(block['gradient_hartree_per_bohr']) * HARTREE_EV / BOHR_A,
                                   rtol=1e-9)
        self.assertLess(np.abs(block['net_force_ev_per_angstrom']).max(), 1e-3)
        s = 1e-4
        volume = block['cell_volume_bohr3']
        numeric_stress = (self.run_cfg(self.strained(cfg, s), 'sp')['total_energy_hartree']
                          - self.run_cfg(self.strained(cfg, -s), 'sm')['total_energy_hartree']) / (2 * s * volume)
        stress = np.asarray(block['stress_hartree_per_bohr3'])
        self.assertAlmostEqual(stress[0, 0], numeric_stress, delta=1e-9)
        np.testing.assert_allclose(block['stress_gpa'], stress * 29421.015697, rtol=1e-6)
        self.assertAlmostEqual(block['pressure_gpa'], -np.trace(block['stress_gpa']) / 3, places=10)

    def test_smearing_forces_are_free_energy_derivatives(self):
        cfg = silicon(calculation='gradient', smearing_ev=0.1)
        block = self.run_cfg(cfg, 'ref')['gradient']
        self.assertEqual(block['derivative_of'], 'Mermin free energy F = E - TS')
        h = 1e-3
        plus, minus = self.run_cfg(self.displaced(cfg, h), 'p'), self.run_cfg(self.displaced(cfg, -h), 'm')
        free = (plus['free_energy_hartree'] - minus['free_energy_hartree']) / (2 * h)
        self.assertAlmostEqual(block['gradient_hartree_per_bohr'][1][0], free, delta=1e-7)

    def test_gdf_and_fftdf_energies_are_never_mixed(self):
        from energy_analysis import reaction_energy
        self.run_cfg(silicon(), 'fft')
        self.run_cfg(silicon(integral_scheme='gdf', fft_mesh=()), 'gdf')
        with self.assertRaisesRegex(ValueError, 'integral_scheme'):
            reaction_energy([(1, self.root / 'fft'), (-1, self.root / 'gdf')])


if __name__ == '__main__':
    unittest.main()
