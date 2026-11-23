import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from engines.frequencies import BOHR_ANGSTROM
from molecular_isotopes import DEUTERIUM_MASS_AMU, WATER_EXAMPLE, compare_isotopes, run_isotopes


def water():
    return json.loads(WATER_EXAMPLE.read_bytes())


def spring_model(elements, coords, masses, hessian):
    return {'calculation': 'frequencies', 'converged': True,
            'config': {'unit': 'bohr', 'atoms': [dict(element=e, position=p) for e, p in zip(elements, coords)]},
            'frequencies': {'masses_amu': masses, 'hessian_hartree_per_bohr2': hessian.tolist()}}


class IsotopeTests(unittest.TestCase):
    def test_h2_hd_d2_reduced_mass_and_imaginary_mode(self):
        stretch = np.array([0, 0, 1, 0, 0, -1])
        hessian = np.outer(stretch, stretch) * 0.4
        h, d = 1.007825, DEUTERIUM_MASS_AMU
        source = spring_model(['H', 'H'], [[0, 0, 0], [0, 0, 1.4]], [h, h], hessian)
        for m1, kwargs in ((d, {'deuterate': True}), (h, {'assignments': [f'2={d}']})):
            with self.subTest(m1=m1):
                result = compare_isotopes(source, **kwargs)
                mode, = result['modes']
                ratio = np.sqrt((h/2) / (m1*d/(m1+d)))
                self.assertAlmostEqual(mode['substituted_cm1']/mode['reference_cm1'], ratio, places=12)
                self.assertAlmostEqual(result['substituted_zpe_kj_mol']/result['reference_zpe_kj_mol'], ratio, places=12)
        source['frequencies']['hessian_hartree_per_bohr2'] = (-hessian).tolist()
        result = compare_isotopes(source, deuterate=True)
        self.assertLess(result['modes'][0]['substituted_cm1'], 0)
        self.assertTrue(any('허수' in msg for msg in result['warnings']))
        self.assertEqual(result['reference_zpe_kj_mol'], 0)
        self.assertEqual(result['substituted_zpe_kj_mol'], 0)

    def test_mode_crossing_follows_displacement(self):
        # Orthogonal stretches of a bent triatomic. The second bond drops below the first.
        bonds = np.array([[-1, 0, 0, 1, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 1, 0]])
        hessian = bonds.T @ np.diag([1, 2]) @ bonds
        source = spring_model(['O', 'H', 'H'], [[0, 0, 0], [1, 0, 0], [0, 1, 0]], [16, 1, 1], hessian)
        source['frequencies']['infrared'] = {'dipole_derivative_e': np.outer([1, 0, 0], bonds[1]).tolist()}
        rows = compare_isotopes(source, assignments=['3=4'])['modes']
        self.assertEqual([row['substituted_mode'] for row in rows], [1, 3, 2])
        self.assertAlmostEqual(rows[1]['shift_cm1'], 0, places=8)
        self.assertLess(rows[2]['shift_cm1'], 0)
        for label in ('reference', 'substituted'):
            self.assertLess(rows[1][f'{label}_ir_km_mol'], 1e-15)
            self.assertGreater(rows[2][f'{label}_ir_km_mol'], 100)

    def test_old_water_input_and_output_preservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'results.json'
            source.write_bytes(WATER_EXAMPLE.read_bytes())
            original = source.read_bytes()
            result, first = run_isotopes(source.parent, tmp, deuterate=True)
            hdo, second = run_isotopes(source, tmp, assignments=[f'2={DEUTERIUM_MASS_AMU}'])
            self.assertNotEqual(first, second)
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(json.loads((first/'isotopes.json').read_bytes()), result)
            with (first/'modes.csv').open(encoding='utf-8', newline='') as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 3)
            modes = result['modes']
            np.testing.assert_allclose([r['reference_cm1'] for r in modes],
                                       water()['frequencies']['frequencies_cm1'], rtol=1e-8)
            self.assertTrue(all(row['shift_cm1'] < 0 for row in modes))
            self.assertLess(result['substituted_zpe_kj_mol'], hdo['substituted_zpe_kj_mol'])
            self.assertLess(hdo['substituted_zpe_kj_mol'], result['reference_zpe_kj_mol'])
            self.assertEqual(hdo['atoms'][2]['substituted_mass_amu'], hdo['atoms'][2]['reference_mass_amu'])
            other_units = water()
            other_units['config']['unit'] = 'bohr'
            for atom in other_units['config']['atoms']:
                atom['position'] = (np.asarray(atom['position']) / BOHR_ANGSTROM).tolist()
            self.assertEqual(modes, compare_isotopes(other_units, deuterate=True)['modes'])

    def test_invalid_inputs(self):
        for tokens in ([], ['4=2'], ['0=2'], ['2=nan'], ['2=0'], ['2=-1'], ['2=2', '2=3'], ['D']):
            with self.subTest(tokens=tokens), self.assertRaises(ValueError):
                compare_isotopes(water(), assignments=tokens)
        for source in (dict(water(), converged=False), dict(water(), calculation='gradient'), dict(water(), system='periodic'),
                       dict(water(), frequencies={}), dict(water(), config={})):
            with self.assertRaises(ValueError):
                compare_isotopes(source, deuterate=True)
        source = water()
        source['frequencies']['reference_max_gradient_hartree_per_bohr'] = 0.1
        self.assertTrue(any('평형 구조' in msg for msg in compare_isotopes(source, deuterate=True)['warnings']))
        source['frequencies']['hessian_hartree_per_bohr2'][0][1] += 1
        with self.assertRaisesRegex(ValueError, '헤시안'):
            compare_isotopes(source, deuterate=True)
        source = water()
        for atom in source['config']['atoms']:
            atom['element'] = 'O'
        with self.assertRaisesRegex(ValueError, '수소'):
            compare_isotopes(source, deuterate=True)

    @unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires PySCF')
    def test_water_isotopologues_match_pyscf(self):
        from pyscf import gto
        from pyscf.hessian import thermo
        source = water()
        atoms = source['config']['atoms']
        mol = gto.M(atom=[(a['element'], a['position']) for a in atoms], unit=source['config']['unit'],
                    basis=source['config']['basis'], verbose=0)
        flat = np.asarray(source['frequencies']['hessian_hartree_per_bohr2'])
        hessian = flat.reshape(mol.natm, 3, mol.natm, 3).transpose(0, 2, 1, 3)
        for kwargs in ({'deuterate': True}, {'assignments': [f'2={DEUTERIUM_MASS_AMU}']}):
            result = compare_isotopes(source, **kwargs)
            for label in ('reference', 'substituted'):
                masses = np.array([atom[f'{label}_mass_amu'] for atom in result['atoms']])
                expected = thermo.harmonic_analysis(mol, hessian, mass=masses, imaginary_freq=False)['freq_wavenumber']
                actual = sorted(row[f'{label}_cm1'] for row in result['modes'])
                np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-3)


if __name__ == '__main__':
    unittest.main()
