import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import trapezoid

from engines.frequencies import harmonic_analysis
from engines.infrared import KM_MOL_PER_E2_AMU, broaden_spectrum, mode_intensities
from molecular_isotopes import WATER_EXAMPLE, WATER_IR_EXAMPLE, compare_isotopes, run_isotopes
from molecular_plotting import analysis_figure


class InfraredTests(unittest.TestCase):
    def test_stretch_uses_inverse_reduced_mass(self):
        stretch = np.array([0., 0., -1., 0., 0., 1.])
        hessian = np.outer(stretch, stretch) / 3
        coords = np.array([[0, 0, -1.], [0, 0, 1.]])
        tensor = np.array([stretch, 2 * stretch, 0 * stretch])
        for masses in ([1., 1.], [1., 2.], [2., 2.]):
            modes = harmonic_analysis(hessian, coords, np.array(masses))
            intensity, = mode_intensities(tensor, modes)
            self.assertAlmostEqual(intensity, 5 * (1/masses[0] + 1/masses[1]) * KM_MOL_PER_E2_AMU, places=9)
        saddle = harmonic_analysis(-hessian, coords, np.ones(2))
        self.assertEqual(mode_intensities(tensor, saddle), [None])
        for invalid in (tensor.T, np.full((3, 6), np.nan), {}, None):
            with self.assertRaises(ValueError):
                mode_intensities(invalid, saddle)

    def test_gaussian_area_width_and_imaginary_modes(self):
        series = [[(1000., 27.), (-50., None)], [(2500., 12.)]]
        for fwhm in (.15, 30.):
            grid, curves = broaden_spectrum(series, fwhm)
            self.assertTrue(np.all(np.diff(grid) > 0))
            np.testing.assert_allclose(trapezoid(curves, x=grid, axis=1), [27, 12], rtol=3e-4)
            half = np.interp(1000 + fwhm / 2, grid, curves[0])
            self.assertAlmostEqual(half / curves[0].max(), 0.5, delta=.001)
        grid, curves = broaden_spectrum([[(-10., None)]], 20.)
        self.assertEqual((grid.size, curves.shape), (0, (1, 0)))
        for invalid in (-2, 0., float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                broaden_spectrum(series, invalid)

    def test_reference_export_and_legacy_input(self):
        sample = json.loads(WATER_IR_EXAMPLE.read_bytes())
        stored = sample['frequencies']['infrared']['intensities_km_mol']
        with tempfile.TemporaryDirectory() as temp:
            result, out = run_isotopes(WATER_IR_EXAMPLE, temp, deuterate=True, fwhm=22)
            self.assertEqual(result['ir_labels'], ['H2O', 'D2O'])
            self.assertEqual(result['ir_fwhm_cm1'], 22)
            np.testing.assert_allclose([r['reference_ir_km_mol'] for r in result['modes']], stored, atol=1e-8)
            grid, ref_curve, sub_curve = np.loadtxt(out/'infrared.csv', delimiter=',', skiprows=1, unpack=True)
            self.assertGreater(len(grid), 100)
            for curve, key in ((ref_curve, 'reference_ir_km_mol'), (sub_curve, 'substituted_ir_km_mol')):
                self.assertAlmostEqual(trapezoid(curve, x=grid), sum(r[key] for r in result['modes']), delta=.03)
            self.assertTrue((out/'infrared.png').read_bytes().startswith(b'\x89PNG'))
            self.assertEqual(json.loads((out/'isotopes.json').read_bytes())['ir_labels'], result['ir_labels'])
            old, old_dir = run_isotopes(WATER_EXAMPLE, temp, deuterate=True)
            self.assertIsNone(old['ir_labels'])
            self.assertFalse((old_dir/'infrared.png').exists())
            self.assertTrue(any('쌍극자 미분' in w for w in old['warnings']))
        sample['frequencies']['infrared'] = {}
        with self.assertRaisesRegex(ValueError, '쌍극자'):
            compare_isotopes(sample, deuterate=True)

    def test_gas_summary_retains_thermo_plot(self):
        sample = json.loads(WATER_IR_EXAMPLE.read_bytes())
        sample.update(calculation='gas_thermo', gas_thermo={'rows': []})
        self.assertEqual(len(analysis_figure(sample).axes), 3)


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires PySCF')
class DipoleTests(unittest.TestCase):
    def test_ionic_translation_and_inactive_stretch(self):
        from molecular import DFTConfig, parse_coordinates
        from engines import calculate
        cfg = DFTConfig(atoms=parse_coordinates('H 0 0 -.5; H 0 0 .5'), calculation='frequencies',
                        engine='native', charge=1, spin=1, xc='pbe', basis='sto-3g', grid_level=2,
                        cube_points=12, threads=1)
        with tempfile.TemporaryDirectory() as temp:
            summary = calculate(cfg, Path(temp), lambda _: None)
        ir = summary['frequencies']['infrared']
        # Summing both atomic derivatives translates the entire ion: dμ_i/dt_j = +1 δ_ij.
        tensor = np.asarray(ir['dipole_derivative_e']).reshape(3, 2, 3)
        np.testing.assert_allclose(tensor.sum(axis=1), np.eye(3), atol=1e-5)
        self.assertLess(ir['intensities_km_mol'][0], 1e-8)

    def test_translated_ecp_and_x2c_ions(self):
        from molecular import DFTConfig, parse_coordinates
        from engines.pyscf_engine import calculate_molecule
        for atom, spin, basis, rel in [('Rb', 0, 'def2-svp', 'none'), ('He', 1, 'sto-3g', 'x2c')]:
            with self.subTest(atom=atom), tempfile.TemporaryDirectory() as temp:
                cfg = DFTConfig(atoms=parse_coordinates(f'{atom} .2 -.3 .4'), unit='bohr', charge=1, spin=spin,
                                basis=basis, relativity=rel, grid_level=2, cube_points=12, threads=1)
                result = calculate_molecule(cfg, Path(temp), lambda _: None, with_dipole=True)
                np.testing.assert_allclose(result['dipole_e_bohr'], [.2, -.3, .4], atol=1e-7)
