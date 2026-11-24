import csv
import json
import math
from pathlib import Path
import tempfile
import unittest

from scipy.constants import R

from gas_mixture import DEFAULT_DATA, DEFAULT_MIX, TEMPERATURES, mixture_report, run_mixture


# NIST-JANAF, tables N-023, O-029 and C-095. Molar Cp in J/(mol K).
JANAF_T = (298.15, 500, 700, 1000, 1200, 2000)
JANAF_CP = {
    'N2': (29.124, 29.580, 30.754, 32.697, 33.723, 35.971),
    'O2': (29.376, 31.091, 32.981, 34.870, 35.667, 37.741),
    'CO2': (37.129, 44.627, 49.564, 54.308, 56.342, 60.350),
}


class MixtureReferenceTests(unittest.TestCase):
    def test_pure_reference_at_tabulated_janaf_points(self):
        for name, tabulated in JANAF_CP.items():
            rows = mixture_report(name + '=1', temperatures=JANAF_T)['reference_comparison']['rows']
            for row, expected in zip(rows, tabulated, strict=True):
                with self.subTest(gas=name, T=row['temperature_K']):
                    # The published Shomate fits approximate the rounded JANAF table.
                    self.assertAlmostEqual(row['nist_cp_J_mol_K'], expected, delta=.032)

    def test_shomate_switches_use_upper_interval_at_the_boundary(self):
        # Left- and right-hand limits evaluated from the WebBook A–E coefficients.
        limits = [('N2', 500, 29.5814235, 29.58134825), ('O2', 700, 32.97548930, 32.97616924),
                  ('CO2', 1200, 56.34605444, 56.31141827)]
        for name, switch, left, right in limits:
            with self.subTest(gas=name):
                result = mixture_report(name + '=1', temperatures=[switch - 1e-6, switch, switch + 1e-6])
                comparison = result['reference_comparison']
                self.assertEqual(comparison['species_sources'][name]['shomate_switch_K'], switch)
                for row, expected in zip(comparison['rows'], (left, right, right), strict=True):
                    self.assertAlmostEqual(row['nist_cp_J_mol_K'], expected, delta=3e-8)

    def test_binary_reference_from_tabulated_cp_and_sound_speed_relation(self):
        mass = (.79 * 28.014 + .21 * 31.998) / 1000
        result = mixture_report('O2=.21,CO2=0,N2=.79', temperatures=JANAF_T)
        comparison = result['reference_comparison']
        self.assertAlmostEqual(result['molar_mass_kg_mol'], mass)
        self.assertEqual(set(comparison['species_sources']), {'N2', 'O2'})
        for row, n2, o2 in zip(comparison['rows'], JANAF_CP['N2'], JANAF_CP['O2'], strict=True):
            with self.subTest(T=row['temperature_K']):
                cp = .79 * n2 + .21 * o2
                speed = math.sqrt(R * row['temperature_K'] / (mass * (1 - R / cp)))
                self.assertAlmostEqual(row['nist_cp_J_mol_K'], cp, delta=.007)
                self.assertAlmostEqual(row['nist_cp_J_kg_K'], cp / mass, delta=.25)
                self.assertAlmostEqual(row['nist_derived_sound_speed_m_s'], speed, delta=.03)
                self.assertEqual((row['x_N2'], row['x_O2'], row['x_CO2'], row['pressure_Pa']), (.79, .21, 0., 100_000))
        ternary = mixture_report('N2=.7,O2=.2,CO2=.1', temperatures=[1200])['reference_comparison']
        self.assertEqual(set(ternary['species_sources']), set(JANAF_CP))
        self.assertAlmostEqual(ternary['rows'][0]['nist_cp_J_mol_K'], .7*33.723 + .2*35.667 + .1*56.342, delta=.007)

    def test_reference_is_independent_of_the_input_nasa7(self):
        data = json.loads(DEFAULT_DATA.read_text(encoding='utf-8'))
        nitrogen = next(item for item in data['species'] if item['name'] == 'N2')
        for block in nitrogen['thermo']['data']:
            block[0] += .75
        temperatures = (298.15, 1000, 1000.01, 2000)
        baseline = mixture_report(DEFAULT_MIX, temperatures=temperatures)['reference_comparison']['rows']
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'custom.yaml'
            path.write_text(json.dumps(data), encoding='utf-8')
            altered = mixture_report(DEFAULT_MIX, files=[path], temperatures=temperatures)['reference_comparison']['rows']
        for old, new in zip(baseline, altered, strict=True):
            for key in ('nist_cp_J_mol_K', 'nist_cp_J_kg_K', 'nist_derived_gamma', 'nist_derived_sound_speed_m_s'):
                self.assertEqual(old[key], new[key])
            self.assertAlmostEqual(new['cp_delta_J_mol_K'] - old['cp_delta_J_mol_K'], .79 * .75 * R)
            self.assertGreater(new['cp_delta_percent'], 0)
            self.assertLess(new['sound_speed_delta_m_s'], 0)
            self.assertAlmostEqual(new['cp_delta_percent'], 100 * (new['nasa7_cp_J_mol_K'] / new['nist_cp_J_mol_K'] - 1))
            self.assertAlmostEqual(new['sound_speed_delta_percent'],
                                   100 * (new['nasa7_sound_speed_m_s'] / new['nist_derived_sound_speed_m_s'] - 1))

    def test_extrema_follow_only_the_sampled_temperatures(self):
        comparison = mixture_report('CO2=1', temperatures=[2000, 1000, 298.15])['reference_comparison']
        for peak in comparison['extrema']:
            self.assertEqual(peak['delta_peak_temperature_K'], 2000)
            self.assertEqual(peak['percent_peak_temperature_K'], 298.15)
            self.assertEqual((peak['x_CO2'], peak['pressure_Pa'], peak['sampled_points']), (1, 100_000, 3))
        single = mixture_report(DEFAULT_MIX, temperatures=[400])['reference_comparison']
        self.assertEqual(single['sampled_temperature_range_K'], [400, 400])
        self.assertEqual(single['valid_temperature_range_K'], [298.15, 2000])
        self.assertEqual(single['sampled_points'], 1)
        self.assertTrue(all(p['delta_peak_temperature_K'] == 400 for p in single['extrema']))
        with self.assertRaises(ValueError):
            mixture_report(DEFAULT_MIX, temperatures=[])

    def test_binary_exports_preserve_previous_results_and_reference_units(self):
        with tempfile.TemporaryDirectory() as directory:
            report, path = run_mixture(DEFAULT_MIX, Path(directory))
            stored = json.loads((path / 'mixture.json').read_text(encoding='utf-8'))['reference_comparison']
            self.assertEqual(stored, report['reference_comparison'])
            for filename, key in [('reference_comparison.csv', 'rows'), ('reference_extrema.csv', 'extrema')]:
                with (path / filename).open(encoding='utf-8', newline='') as stream:
                    entries = list(csv.DictReader(stream))
                self.assertEqual(entries, [{k: str(v) for k, v in row.items()} for row in stored[key]])
            self.assertEqual(stored['sampled_points'], len(TEMPERATURES))
            self.assertAlmostEqual(stored['rows'][-1]['cp_delta_J_kg_K'], -16.3539661, places=6)
            self.assertAlmostEqual(stored['rows'][-1]['sound_speed_delta_m_s'], 1.6909527, places=6)
            for png in ('mixture.png', 'reference_comparison.png'):
                self.assertTrue((path / png).read_bytes().startswith(b'\x89PNG\r\n\x1a\n'))
            with DEFAULT_DATA.with_name('air_sound_speed.csv').open(encoding='utf-8', newline='') as stream:
                previous = list(csv.DictReader(stream))
            for row, old in zip(report['rows'], previous, strict=True):
                self.assertEqual(row.keys(), old.keys())
                for key, value in row.items():
                    self.assertAlmostEqual(value, float(old[key]), places=9)
