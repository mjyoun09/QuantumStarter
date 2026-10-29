"""Climbing-image NEB: input rules, symmetric NH3 inversion path, agreement with the TS search."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from molecular import Atom, DFTConfig, coordinates_text, parse_coordinates
from scf_control import signature

NH3 = 'N 0 0 0.12; H 0 0.94 -0.27; H 0.814 -0.47 -0.27; H -0.814 -0.47 -0.27'
READY = bool(importlib.util.find_spec('pyscf') and importlib.util.find_spec('geometric'))


def mirrored(atoms):
    return tuple(Atom(a.element, (a.position[0], a.position[1], -a.position[2])) for a in atoms)


class NebInputTests(unittest.TestCase):
    def test_product_rules(self):
        atoms = parse_coordinates(NH3)
        cfg = DFTConfig(atoms=atoms, calculation='neb', neb_product=coordinates_text(mirrored(atoms)))
        self.assertEqual(DFTConfig.from_dict(cfg.to_dict()), cfg)
        self.assertEqual(signature(cfg), signature(replace(cfg, calculation='single_point', neb_product='')))
        for kwargs in ({'neb_product': ''}, {'neb_product': 'H 0 0 0; N 0 0 1; H 1 0 0; H 0 1 0'},
                       {'neb_product': 'N 0 0 0; H 0 0 1; H 1 0 0'}, {'neb_images': 2}, {'neb_images': 42},
                       {'neb_maxg_ev_per_angstrom': 0}):
            with self.subTest(**kwargs), self.assertRaises(ValueError):
                replace(cfg, **kwargs)
        with self.assertRaises(ValueError):
            replace(cfg, calculation='optimize')  # product only for NEB


@unittest.skipUnless(READY, 'Requires WSL PySCF and geomeTRIC')
class NebPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.base = DFTConfig(atoms=parse_coordinates(NH3), xc='pbe', basis='sto-3g', engine='native',
                              grid_level=3, cube_points=12, threads=1)

    def tearDown(self):
        self.temp.cleanup()

    def test_ammonia_inversion_path_and_ts_agreement(self):
        from engines import calculate
        minimum = calculate(replace(self.base, calculation='optimize', opt_convergence='tight'), self.root / 'min',
                            lambda _: None)
        reactant = DFTConfig.from_dict(minimum['final_geometry_config']).atoms
        cfg = replace(self.base, atoms=reactant, calculation='neb', neb_product=coordinates_text(mirrored(reactant)),
                      neb_images=7)
        result = calculate(cfg, self.root / 'neb', lambda _: None)
        block = result['neb']
        self.assertEqual(result['status'], 'converged')
        self.assertTrue(block['converged'])
        self.assertLess(block['max_gradient_ev_per_angstrom'], 0.05)
        profile = [row['relative_kj_mol'] for row in block['image_data']]
        np.testing.assert_allclose(profile, profile[::-1], atol=1e-3)  # mirror-symmetric path
        self.assertEqual(block['highest_image'], 3)
        top = np.array([a['position'] for a in result['final_geometry_config']['atoms']])
        self.assertLess(abs(top[0, 2] - top[1:, 2].mean()), 1e-3)  # planar
        refined = calculate(replace(self.base, atoms=DFTConfig.from_dict(result['final_geometry_config']).atoms,
                                    calculation='transition_state'), self.root / 'ts', lambda _: None)
        barrier = (refined['total_energy_hartree'] - minimum['total_energy_hartree']) * 2625.4996394799
        self.assertAlmostEqual(block['barrier_from_reactant_kj_mol'], barrier, delta=0.1)
        run = self.root / 'neb'
        self.assertEqual(len((run / 'neb.csv').read_text(encoding='utf-8').splitlines()), 8)
        self.assertEqual((run / 'neb_path.xyz').read_text(encoding='utf-8').count('image '), 7)
        self.assertEqual(len(list((run / 'neb').iterdir())), block['gradient_calculations'])

    def test_missing_library_is_explicit(self):
        from engines import calculate
        cfg = replace(self.base, calculation='neb', neb_product=coordinates_text(mirrored(self.base.atoms)))
        with patch.dict(sys.modules, {'geometric': None}):
            with self.assertRaisesRegex(RuntimeError, 'geomeTRIC'):
                calculate(cfg, self.root / 'missing', lambda _: None)


if __name__ == '__main__':
    unittest.main()
