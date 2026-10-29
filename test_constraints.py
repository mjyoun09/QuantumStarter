"""Constrained optimization, relaxed scans and transition-state search (geomeTRIC + native forces)."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from molecular import DFTConfig, parse_coordinates, parse_scan
from scf_control import signature

WATER = 'O 0 0 0; H 0 -.80 .55; H 0 .80 .55'
READY = bool(importlib.util.find_spec('pyscf') and importlib.util.find_spec('geometric'))


class ConstraintInputTests(unittest.TestCase):
    def test_scan_spec_and_field_rules(self):
        self.assertEqual(parse_scan('distance 1 2 0.6 1.0 5', 2), ('distance', (1, 2), [0.6, 0.7, 0.8, 0.9, 1.0]))
        self.assertEqual(parse_scan('angle 2 1 3 90 120 4', 3)[2], [90.0, 100.0, 110.0, 120.0])
        for spec in ('distance 1 2 0.6 1.0', 'distance 1 1 0.6 1.0 5', 'distance 1 3 0.6 1.0 5',
                     'torsion 1 2 3 4 0 90 4', 'distance 1 2 0.6 0.6 5', 'distance 1 2 0.6 1.0 1'):
            with self.subTest(spec=spec), self.assertRaises(ValueError):
                parse_scan(spec, 2)
        cfg = DFTConfig(calculation='optimize', opt_constraints='$freeze;xyz 1')
        self.assertEqual(signature(cfg), signature(DFTConfig()))
        with self.assertRaises(ValueError):
            DFTConfig(opt_constraints='$freeze;xyz 1')  # single_point has no geometry to constrain
        with self.assertRaises(ValueError):
            DFTConfig(calculation='optimize', scan_spec='distance 1 2 0.6 1.0 5')
        DFTConfig(calculation='scan', scan_spec='angle 2 1 3 90 120 4')


@unittest.skipUnless(READY, 'Requires WSL PySCF and geomeTRIC')
class ConstraintPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.cfg = DFTConfig(atoms=parse_coordinates(WATER), xc='pbe', basis='sto-3g', engine='native',
                             grid_level=3, cube_points=12, threads=1, calculation='optimize', opt_convergence='tight')

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def final(result):
        return np.array([a['position'] for a in result['final_geometry_config']['atoms']])

    def test_angle_set_and_atom_freeze(self):
        from engines import calculate
        result = calculate(replace(self.cfg, opt_constraints='$set;angle 2 1 3 100'), self.root / 'angle', lambda _: None)
        x = self.final(result)
        a, b = x[1] - x[0], x[2] - x[0]
        self.assertEqual(result['status'], 'converged')
        self.assertAlmostEqual(np.degrees(np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))), 100, delta=1e-6)
        self.assertAlmostEqual(np.linalg.norm(a), np.linalg.norm(b), delta=1e-5)
        self.assertIsNone(result['optimization']['final_gradient_meets_force_criteria'])
        self.assertTrue((self.root / 'angle' / 'constraints.txt').is_file())
        frozen = calculate(replace(self.cfg, opt_constraints='$freeze;xyz 1'), self.root / 'freeze', lambda _: None)
        # geomeTRIC enforces constraints to its convergence tolerance, not exactly.
        self.assertLess(np.linalg.norm(self.final(frozen)[0] - np.zeros(3)), 1e-3)
        with self.assertRaisesRegex(ValueError, '제약'):
            calculate(replace(self.cfg, opt_constraints='$set;angle 2 1 100'), self.root / 'bad', lambda _: None)

    def test_relaxed_scan_matches_single_points(self):
        from engines import calculate
        scan = replace(self.cfg, atoms=parse_coordinates('H 0 0 0; H 0 0 .74'), calculation='scan',
                       scan_spec='distance 1 2 0.6 1.0 5')
        result = calculate(scan, self.root / 'scan', lambda _: None)
        self.assertEqual(result['status'], 'converged')
        rows = result['scan']['points']
        self.assertEqual(len(rows), 5)
        for row in rows:
            c = np.array(row['coordinates_angstrom'])
            self.assertAlmostEqual(np.linalg.norm(c[1] - c[0]), row['value'], delta=1e-6)
            single = calculate(replace(scan, calculation='single_point', scan_spec='',
                                       atoms=parse_coordinates(f'H 0 0 0; H 0 0 {row["value"]}')),
                               self.root / f'single{row["point"]}', lambda _: None)
            self.assertAlmostEqual(row['energy_hartree'], single['total_energy_hartree'], delta=1e-8)
        self.assertEqual(result['scan']['minimum_point'], 2)
        self.assertEqual(len((self.root / 'scan' / 'scan.csv').read_text(encoding='utf-8').splitlines()), 6)
        self.assertTrue((self.root / 'scan' / 'density.cube').is_file())

    def test_ammonia_inversion_transition_state(self):
        from engines import calculate
        ts = replace(self.cfg, atoms=parse_coordinates('N 0 0 0.05; H 0 1.01 -0.01; H .8747 -.505 -0.01; '
                                                       'H -.8747 -.505 -0.01'),
                     calculation='transition_state', opt_convergence='default')
        result = calculate(ts, self.root / 'ts', lambda _: None)
        self.assertEqual(result['status'], 'converged')
        self.assertEqual(result['optimization']['hessian_gradient_calculations'], 24)
        x = self.final(result)
        self.assertLess(abs(x[0, 2] - x[1:, 2].mean()), 2e-3)  # planar
        frequencies = calculate(replace(ts, atoms=tuple(DFTConfig.from_dict(result['final_geometry_config']).atoms),
                                        calculation='frequencies'), self.root / 'freq', lambda _: None)['frequencies']
        self.assertEqual(frequencies['stationary_point_character'], 'transition_state_candidate')
        self.assertEqual(len(frequencies['imaginary_frequencies_cm1']), 1)
        mode = np.asarray(frequencies['normal_modes_cartesian'][0])
        self.assertGreater(np.linalg.norm(mode[:, 2]), .99)  # umbrella motion out of the plane


if __name__ == '__main__':
    unittest.main()
