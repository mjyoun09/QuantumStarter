import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from dft_workflow import CalculationError, RunControl
from gas_mixture import DEFAULT_MIX, mixture_report, parse_fractions
from molecular import DFTConfig
from molecular_convergence import difference, run_comparison


class MixtureTests(unittest.TestCase):
    def test_composition_and_temperature_bounds(self):
        self.assertEqual(parse_fractions(' n2=1, O2=0 '), {'N2': 1.})
        for mix in ('N2=79,O2=21', 'N2=.5', 'N2=.5,n2=.5', 'CO=1', 'O2=nan', 'CO2=-1,N2=2'):
            with self.subTest(mix=mix), self.assertRaises(ValueError):
                parse_fractions(mix)
        for t in (298., 2000.1, float('nan')):
            with self.assertRaises(ValueError):
                mixture_report(DEFAULT_MIX, temperatures=[t])

    def test_molar_mixing_and_mass_conversion(self):
        def row(mix):
            return mixture_report(mix, temperatures=[1200])['rows'][0]
        mix = row(DEFAULT_MIX)
        for prop in ('cp_J_mol_K','cv_J_mol_K'):
            self.assertAlmostEqual(mix[prop], .79*row('N2=1')[prop]+.21*row('O2=1')[prop])
        self.assertAlmostEqual(mix['cp_J_kg_K']*(.79*.028014+.21*.031998), mix['cp_J_mol_K'])


class ComparisonTests(unittest.TestCase):
    def test_changed_crossing_count(self):
        a = dict(basis='a', grid=1, values={('crossings',0,'count'):1, ('crossing_midpoint',1,'angstrom'):1.5})
        b = dict(basis='b', grid=1, values={('crossings',0,'count'):2, ('crossing_midpoint',1,'angstrom'):1.3,
                 ('crossing_midpoint',2,'angstrom'):2.1})
        rows = {r['observable']+str(r['rank']):r for r in difference('basis', a, b)}
        self.assertEqual(rows['crossings0']['delta'], 1)
        self.assertEqual(rows['crossing_midpoint1']['status'], 'different_crossing_count')
        self.assertIsNone(rows['crossing_midpoint2']['delta'])

    def test_preserve_partial_and_cancelled_runs(self):
        cfg = DFTConfig(engine='native', calculation='uvvis')
        def run(c, root, backend, control):
            path = root/f'{c.basis}_{c.grid_level}'
            if c.basis == 'sto-3g' and c.grid_level == 2:
                raise CalculationError('failed SCF', path)
            return {'status':'converged', 'uvvis':{'peaks':[dict(state=1, energy_eV=c.grid_level, oscillator_strength=0.)]}}, path
        with tempfile.TemporaryDirectory() as temp, patch('molecular_convergence.run_calculation', side_effect=run):
            opts = dict(config=cfg, bases=['sto-3g','def2-svp'], grids=[1,2], output_root=Path(temp), progress=lambda _: None)
            report, path = run_comparison(**opts)
            self.assertEqual(report['status'], 'incomplete')
            self.assertEqual(len(report['runs']), 4)
            self.assertEqual(len(report['differences']), 4)
            self.assertEqual(json.loads((path/'comparison.json').read_text(encoding='utf-8'))['runs'][1]['status'], 'error')
            control = RunControl()
            control.cancel()
            cancelled, second = run_comparison(**opts, control=control)
            self.assertEqual(cancelled['status'], 'cancelled')
            self.assertFalse(cancelled['runs'])
            self.assertNotEqual(second, path)
