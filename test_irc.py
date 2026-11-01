"""IRC from a transition state: input rules, best-fit helpers, NH3 inversion both ways, rejection of a minimum."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from molecular import DFTConfig, parse_coordinates
from scf_control import signature

# NH3 inversion TS from this repository's transition_state search at the same level
# (native PBE/STO-3G, grid 3, default criteria), rounded to 1e-6 Å. Max force 2.9e-3 eV/Å.
NH3_TS = ('N -0.001007 0.269182 0.016025; H -0.001629 1.295287 -0.003789; '
          'H 0.887857 -0.243738 0.026694; H -0.889597 -0.244210 0.024835')
NH3 = 'N 0 0 0.12; H 0 0.94 -0.27; H 0.814 -0.47 -0.27; H -0.814 -0.47 -0.27'
READY = bool(importlib.util.find_spec('pyscf') and importlib.util.find_spec('geometric'))


def rotation(axis, angle):
    axis = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


class IrcInputTests(unittest.TestCase):
    def test_fields(self):
        cfg = DFTConfig(atoms=parse_coordinates(NH3_TS), calculation='irc')
        self.assertEqual(signature(cfg), signature(replace(cfg, calculation='single_point')))
        with self.assertRaisesRegex(ValueError, 'IRC'):
            DFTConfig(atoms=parse_coordinates('H 0 0 0'), calculation='irc')
        with self.assertRaisesRegex(ValueError, '제약'):
            DFTConfig(atoms=parse_coordinates(NH3_TS), calculation='irc', opt_constraints='$freeze;xyz 1')

    def test_best_fit_distances(self):
        from engines.irc import rmsd, superposition_distance
        x = np.array([a.position for a in parse_coordinates(NH3)])
        moved = x @ rotation([1, 2, 3], 0.7).T + [0.3, -1.2, 2.0]
        self.assertLess(rmsd(x, moved), 1e-12)
        mirror = x * [1, 1, -1]  # Labelled mirror image: no proper rotation superposes it.
        self.assertGreater(rmsd(x, mirror), 0.1)
        self.assertLess(rmsd(x, mirror @ rotation([0, 1, 1], 1.1).T, allow_reflection=True), 1e-12)
        # A pure displacement of one atom, after the best fit, is never longer than the displacement.
        masses = np.array([14.007, 1.008, 1.008, 1.008])
        shifted = x.copy()
        shifted[1] += [0, 0.01, 0]
        self.assertLessEqual(superposition_distance(x, shifted, masses), np.sqrt(1.008) * 0.01 + 1e-15)
        self.assertGreater(superposition_distance(x, shifted, masses), 0.5 * np.sqrt(1.008) * 0.01)


@unittest.skipUnless(READY, 'Requires WSL PySCF and geomeTRIC')
class IrcPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.base = DFTConfig(atoms=parse_coordinates(NH3_TS), xc='pbe', basis='sto-3g', engine='native',
                              grid_level=3, cube_points=12, threads=1, calculation='irc')

    def tearDown(self):
        self.temp.cleanup()

    def test_ammonia_inversion_both_directions(self):
        from engines import calculate
        from engines.irc import rmsd
        result = calculate(self.base, self.root / 'irc', lambda _: None)
        block = result['irc']
        self.assertEqual(result['status'], 'converged')
        self.assertTrue(block['converged'])
        self.assertEqual(len([f for f in block['ts_hessian']['imaginary_frequencies_cm1'] if f < -20]), 1)
        self.assertEqual(block['hessian_gradient_calculations'], 24)
        frames, ends = block['frames'], block['endpoints']
        ts = block['ts_frame']
        self.assertEqual(frames[ts]['branch'], 'ts')
        self.assertEqual(frames[ts]['s_bohr_amu05'], 0.0)
        self.assertTrue(np.all(np.diff([f['s_bohr_amu05'] for f in frames]) > 0))
        # Steepest descent: the energy never rises while walking away from the TS.
        self.assertLess(max(block['max_uphill_hartree'].values()), 1e-8)
        for end in ends.values():
            self.assertTrue(end['meets_gradient_criteria'] and end['moved_from_ts'])
            self.assertGreater(end['barrier_to_ts_kj_mol'], 10)
            x = np.array(frames[end['frame']]['coordinates_angstrom'])
            normal = np.cross(x[2] - x[1], x[3] - x[1])
            self.assertGreater(abs((x[0] - x[1]) @ normal / np.linalg.norm(normal)), 0.3)  # pyramidal
        # Both ends are the same minimum, mirror images of each other (labelled atoms).
        self.assertAlmostEqual(ends['forward']['energy_hartree'], ends['backward']['energy_hartree'], delta=2e-6)
        self.assertLess(block['endpoint_rmsd_with_reflection_angstrom'], 2e-3)
        self.assertGreater(block['endpoint_rmsd_angstrom'], 0.1)
        self.assertAlmostEqual(ends['forward']['path_length_bohr_amu05'], ends['backward']['path_length_bohr_amu05'],
                               delta=0.02 * ends['forward']['path_length_bohr_amu05'])
        self.assertTrue(any('거울상' in w for w in result['warnings']))
        # The top level is the TS recomputed with the requested output settings.
        self.assertAlmostEqual(result['total_energy_hartree'], frames[ts]['energy_hartree'], delta=1e-8)
        # An independent minimization reaches the same minimum as the IRC ends.
        minimum = calculate(replace(self.base, atoms=parse_coordinates(NH3), calculation='optimize'),
                            self.root / 'min', lambda _: None)
        self.assertAlmostEqual(minimum['total_energy_hartree'], ends['forward']['energy_hartree'], delta=2e-6)
        optimized = [a['position'] for a in minimum['final_geometry_config']['atoms']]
        self.assertLess(rmsd(optimized, frames[0]['coordinates_angstrom'], allow_reflection=True), 5e-3)
        run = self.root / 'irc'
        self.assertEqual((run / 'irc_path.xyz').read_text(encoding='utf-8').count('frame '), len(frames))
        self.assertEqual(len((run / 'irc.csv').read_text(encoding='utf-8').splitlines()), len(frames) + 1)
        for name in ('irc_forward_end.xyz', 'irc_backward_end.xyz', 'final_settings.json', 'density.cube'):
            self.assertTrue((run / name).is_file(), name)
        self.assertFalse((run / 'geometric.tmp').exists())
        self.assertEqual(len(list((run / 'steps').iterdir())), block['gradient_calculations'])

    def test_step_limit_recovers_path_and_minimizes_ends(self):
        # 12 IRC steps per direction stop short of both minima (the full path needs ~15).
        # geomeTRIC hands the forward direction over silently and raises in the backward one;
        # the path comes back from its xyz file and both ends are finished by minimization.
        from engines import calculate
        result = calculate(replace(self.base, opt_max_steps=12), self.root / 'short', lambda _: None)
        block = result['irc']
        self.assertFalse(block['geometric_converged'])
        self.assertTrue(block['path_source'].startswith('geometric_irc.xyz'))
        self.assertEqual(result['status'], 'converged')
        frames, ends = block['frames'], block['endpoints']
        self.assertTrue(np.all(np.diff([f['s_bohr_amu05'] for f in frames]) > 0))
        self.assertAlmostEqual(frames[block['ts_frame']]['energy_hartree'], result['total_energy_hartree'], delta=1e-8)
        for name, end in ends.items():
            relaxation = end['relaxation']
            self.assertTrue(relaxation['converged'] and end['meets_gradient_criteria'])
            self.assertFalse(end['irc_end']['meets_gradient_criteria'])
            self.assertLess(relaxation['energy_change_hartree'], 0)
            self.assertTrue((self.root / 'short' / f'endpoint_{name}' / 'optimized.xyz').is_file())
        self.assertAlmostEqual(ends['forward']['energy_hartree'], ends['backward']['energy_hartree'], delta=2e-6)
        self.assertLess(block['endpoint_rmsd_with_reflection_angstrom'], 5e-3)
        self.assertTrue(sum('일반 최소화로 끝점을 마무리' in w for w in result['warnings']) == 2)

    def test_minimum_is_not_a_starting_point(self):
        from engines import calculate
        with self.assertRaisesRegex(ValueError, '1차 안장점'):
            calculate(replace(self.base, atoms=parse_coordinates(NH3)), self.root / 'minimum', lambda _: None)

    def test_missing_library_is_explicit(self):
        from engines import calculate
        with patch.dict(sys.modules, {'geometric': None}):
            with self.assertRaisesRegex(RuntimeError, 'geomeTRIC'):
                calculate(self.base, self.root / 'missing', lambda _: None)


if __name__ == '__main__':
    unittest.main()
