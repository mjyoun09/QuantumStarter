"""Periodic DOS and band structure: formula checks on synthetic data, real silicon consistency."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from engines.periodic_analysis import density_of_states, frontier

HARTREE_TO_EV = 27.211386245988


class DosFormulaTests(unittest.TestCase):
    def test_restricted_and_unrestricted_counts(self):
        energies = np.array([[-0.5, -0.2, 0.1, 0.3], [-0.45, -0.25, 0.15, 0.35]])  # 2 k-points, Hartree
        occupations = np.array([[2, 2, 0, 0], [2, 2, 0, 0]], dtype=float)
        result = density_of_states(energies, occupations, [.5, .5], False, 0.05, HARTREE_TO_EV)
        self.assertAlmostEqual(result['integrated_states'], 8.0, places=4)  # 4 bands x 2 spins
        self.assertEqual(result['expected_states'], 8.0)
        self.assertAlmostEqual(result['occupied_electrons_from_weights'], 4.0)
        self.assertAlmostEqual(result['homo_ev'], -0.2 * HARTREE_TO_EV)
        self.assertAlmostEqual(result['lumo_ev'], 0.1 * HARTREE_TO_EV)
        self.assertAlmostEqual(result['mesh_gap_ev'], 0.3 * HARTREE_TO_EV)
        spin = density_of_states(np.stack((energies, energies + .01)), np.stack((occupations / 2, occupations / 2)),
                                 [.5, .5], True, 0.05, HARTREE_TO_EV)
        self.assertEqual(len(spin['dos_by_spin']), 2)
        self.assertAlmostEqual(spin['integrated_states'], 8.0, places=4)

    def test_metal_has_no_gap(self):
        energies = np.array([[-0.2, -0.01, 0.01]])
        occupations = np.array([[2, 1.2, 0.8]])  # smeared: the 1.2 state counts as occupied
        homo, lumo = frontier(energies, occupations, False)
        self.assertEqual((homo, lumo), (-0.01, 0.01))
        self.assertAlmostEqual(density_of_states(energies, occupations, [1.], False, .05, HARTREE_TO_EV)
                               ['occupied_electrons_from_weights'], 4.0)


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class SiliconBandTests(unittest.TestCase):
    def test_gamma_silicon_dos_and_bands(self):
        from engines import calculate
        from periodic import load_calculation_config
        cfg = replace(load_calculation_config('examples/silicon_bulk.json'), kmesh=(1, 1, 1), smearing_ev=0.0,
                      conv_tol=1e-10, band_path='GXL', band_points=12, threads=1)
        with tempfile.TemporaryDirectory() as temp:
            result = calculate(cfg, temp, lambda _: None)
            self.assertTrue((Path(temp) / 'dos.csv').is_file() and (Path(temp) / 'bands.csv').is_file())
        dos, bands = result['dos'], result['band_structure']
        self.assertEqual(result['status'], 'converged')
        self.assertAlmostEqual(dos['integrated_states'], dos['expected_states'], delta=1e-3)
        self.assertAlmostEqual(dos['occupied_electrons_from_weights'], result['explicit_electrons'], places=8)
        self.assertEqual(bands['special_labels'][0], 'G')
        self.assertEqual(len(bands['axis_inverse_angstrom']), 12)
        # Bands at the SCF mesh reproduce the SCF eigenvalues (PySCF get_bands with GDF: ~1 meV).
        self.assertLess(bands['mesh_consistency_max_ev'], 5e-3)
        gap = bands['gaps_by_spin'][0]
        self.assertGreater(gap['direct_gap_ev'], gap['gap_ev'] - 1e-9)


if __name__ == '__main__':
    unittest.main()
