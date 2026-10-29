"""Reaction bookkeeping and refusal of incompatible/unconverged data."""
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from molecular import DFTConfig, parse_coordinates
from energy_analysis import adsorption_energy, reaction_energy, compare_runs, save_analysis


class EnergyAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.hconfig = DFTConfig(atoms=parse_coordinates('H 0 0 0'), spin=1, basis='sto-3g', engine='native')
        self.h = self.fixture('h', self.hconfig, -.5)
        self.h2 = self.fixture('h2', replace(self.hconfig, atoms=parse_coordinates('H 0 0 0; H 0 0 .74'), spin=0), -1.1)

    def tearDown(self):
        self.temp.cleanup()

    def fixture(self, name, config, energy, **extra):
        path = self.root / (name + '.json')
        data = dict(config=config.to_dict(), total_energy_hartree=energy, converged=True,
                    status='converged', engine='quantumstarter-native', formula=config.formula)
        data.update(extra)
        path.write_text(json.dumps(data), encoding='utf-8')
        return path

    def test_sign_stoichiometry_sources_and_export(self):
        report = adsorption_energy(self.h2, self.h, self.h)
        self.assertAlmostEqual(report['energy_hartree'], -.1)
        self.assertEqual(len(report['terms']), 3)
        self.assertEqual(len(report['terms'][0]['source']['sha256']), 64)
        half = reaction_energy([(.5, self.h2), (-1, self.h)])
        self.assertAlmostEqual(half['energy_hartree'], -.05)
        self.assertTrue(save_analysis(report, self.root / 'analysis').is_file())
        self.assertEqual(compare_runs([self.h, self.h2])['analysis'], 'run_comparison')

    def test_incompatible_nonconverged_and_imbalance_rejected(self):
        for field, value in [('basis', 'def2-svp'), ('xc', 'pbe0'), ('grid_level', 4), ('conv_tol', 1e-6)]:
            changed = self.fixture(field, replace(self.hconfig, **{field: value}), -.5)
            with self.assertRaisesRegex(ValueError, '조건 불일치'):
                adsorption_energy(self.h2, self.h, changed)
        bad = self.fixture('bad', self.hconfig, -.5, status='not_converged', converged=False)
        with self.assertRaisesRegex(ValueError, '수렴'):
            adsorption_energy(self.h2, self.h, bad)
        with self.assertRaisesRegex(ValueError, '수지'):
            adsorption_energy(self.h2, self.h, self.h, factor=.5)
        with self.assertRaises(ValueError):
            adsorption_energy(self.h2, self.h, self.h, factor=float('nan'))
        mixed = self.fixture('external', self.hconfig, -.5, engine='pyscf-molecular')
        with self.assertRaisesRegex(ValueError, 'engine'):
            adsorption_energy(self.h2, self.h, mixed)

    def test_per_element_basis_must_match_across_terms(self):
        # 2 H2 + O2 -> 2 H2O: H stays sto-3g everywhere, O uses def2-svp everywhere.
        water = replace(self.hconfig, atoms=parse_coordinates('O 0 0 0; H 0 -.757 .587; H 0 .757 .587'),
                        spin=0, basis_overrides={'O': 'def2-svp'})
        oxygen = replace(self.hconfig, atoms=parse_coordinates('O 0 0 0; O 0 0 1.21'), spin=2, basis='def2-svp')
        terms = [(2, self.fixture('water', water, -75.)), (-2, self.h2), (-1, self.fixture('o2', oxygen, -149.))]
        report = reaction_energy(terms)
        self.assertAlmostEqual(report['energy_hartree'], -150 + 2.2 + 149)
        self.assertEqual(report['conditions']['basis_by_element'], {'H': 'sto-3g', 'O': 'def2-svp'})
        # The same reaction with def2-svp on H in water only is inconsistent for H.
        water_h = self.fixture('water_h', replace(water, basis='def2-svp', basis_overrides={}), -75.)
        with self.assertRaisesRegex(ValueError, '원소 H'):
            reaction_energy([(2, water_h), (-2, self.h2), (-1, terms[2][1])])


if __name__ == '__main__':
    unittest.main()
