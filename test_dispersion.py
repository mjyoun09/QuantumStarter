"""DFT-D4 dispersion: input contract, library contract, physical limits, SCF integration."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from molecular import ATOMIC_NUMBERS, D4_ELEMENTS, D4_FUNCTIONALS, DFTConfig, ELEMENTS, parse_coordinates
from scf_control import signature
from energy_analysis import reaction_energy

WATER = 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587'
HAS_D4 = bool(importlib.util.find_spec('dftd4') and importlib.util.find_spec('pyscf'))


class DispersionInputTests(unittest.TestCase):
    def test_defaults_choices_and_legacy_files(self):
        self.assertEqual(DFTConfig().dispersion, 'none')
        cfg = DFTConfig(dispersion='D4')
        self.assertEqual(cfg.dispersion, 'd4')
        self.assertEqual(DFTConfig.from_dict(json.loads(json.dumps(cfg.to_dict()))), cfg)
        legacy = DFTConfig().to_dict()
        legacy.pop('dispersion')
        self.assertEqual(DFTConfig.from_dict(legacy), DFTConfig())
        with self.assertRaises(ValueError):
            DFTConfig(dispersion='d3')
        # Dispersion does not change the density, so checkpoints stay interchangeable.
        self.assertEqual(signature(cfg), signature(DFTConfig()))
        native = DFTConfig(engine='native', basis='sto-3g', dispersion='d4')
        self.assertEqual(native.dispersion, 'd4')

    def test_functional_and_element_limits(self):
        with self.assertRaisesRegex(ValueError, 'lda'):
            DFTConfig(xc='lda,vwn', dispersion='d4')
        for xc in sorted(D4_FUNCTIONALS):
            DFTConfig(xc=xc, dispersion='d4')
        with self.assertRaisesRegex(ValueError, 'Rf'):
            DFTConfig(atoms=parse_coordinates('Rf 0 0 0'), basis='dyall-v2z', spin=0, dispersion='d4')
        DFTConfig(atoms=parse_coordinates('Og 0 0 0'), basis='dyall-v2z', dispersion='d4')
        self.assertEqual(len(D4_ELEMENTS), 110)
        self.assertEqual(sorted(set(ELEMENTS) - D4_ELEMENTS, key=ATOMIC_NUMBERS.get),
                         ['Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg'])

    def test_reaction_analysis_refuses_mixed_dispersion(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            base = DFTConfig(atoms=parse_coordinates('H 0 0 0'), spin=1, basis='sto-3g', xc='pbe')
            details = {'model': 'D4', 'parameters': {'s6': 1.0}, 'library': 'dftd4', 'library_version': '4.3.1'}

            def fixture(name, config, energy, dispersion=None):
                path = root / (name + '.json')
                path.write_text(json.dumps(dict(config=config.to_dict(), total_energy_hartree=energy, converged=True,
                                                status='converged', engine='pyscf-molecular', formula=config.formula,
                                                dispersion=dispersion)), encoding='utf-8')
                return path

            h2 = replace(base, atoms=parse_coordinates('H 0 0 0; H 0 0 .74'), spin=0)
            plain = [(1, fixture('h2', h2, -1.1)), (-2, fixture('h', base, -.5))]
            self.assertAlmostEqual(reaction_energy(plain)['energy_hartree'], -.1)
            d4 = [(1, fixture('h2d', replace(h2, dispersion='d4'), -1.1, details)),
                  (-2, fixture('hd', replace(base, dispersion='d4'), -.5, details))]
            self.assertEqual(reaction_energy(d4)['conditions']['dispersion'], 'd4')
            with self.assertRaisesRegex(ValueError, 'dispersion'):
                reaction_energy([d4[0], plain[1]])
            other = dict(details, library_version='9.9')
            with self.assertRaisesRegex(ValueError, 'dispersion_library'):
                reaction_energy([d4[0], (-2, fixture('hv', replace(base, dispersion='d4'), -.5, other))])


@unittest.skipUnless(HAS_D4, 'Requires WSL PySCF and dftd4')
class D4PhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_library_contract(self):
        from dftd4.interface import DampingParam, DispersionModel
        from dftd4.parameters import get_damping_param
        supported = set()
        for z, symbol in enumerate(ELEMENTS, 1):
            try:
                DispersionModel(np.array([z, 1]), np.array([[0., 0, 0], [0, 0, 3]])).get_dispersion(
                    DampingParam(method='pbe'), grad=False)
                supported.add(symbol)
            except RuntimeError:
                pass
        self.assertEqual(supported, set(D4_ELEMENTS))
        for xc in ('lda', 'svwn'):
            with self.assertRaises(KeyError):
                get_damping_param(xc)
        # Explicit published parameters reproduce the library's own method lookup.
        from engines.dispersion import d4_energy
        numbers, coords = np.array([8, 1, 1]), np.array([[0., 0, 0], [0, -1.43, 1.11], [0, 1.43, 1.11]])
        for xc in sorted(D4_FUNCTIONALS):
            result = d4_energy(numbers, coords, 0, xc)
            by_name = DispersionModel(numbers, coords).get_dispersion(DampingParam(method=xc), grad=False)['energy']
            self.assertAlmostEqual(result['energy_hartree'], by_name, delta=1e-14)
            self.assertEqual(result['parameters'], {k: float(v) for k, v in get_damping_param(xc).items()})

    def test_two_body_formula_no_truncation_invariance_and_charge(self):
        from dftd4 import data
        from dftd4.interface import DispersionModel
        from dftd4.parameters import get_damping_param
        from engines.dispersion import d4_energy
        # A dimer has no three-body term. Re-evaluate the D4 rational damping by hand:
        # E = -C6 [s6/(R^6+f^6) + s8 rr/(R^8+f^8)], rr = 3 q_A q_B, f = a1 sqrt(rr) + a2.
        # C6 (charge/CN dependent) comes from the library; units, damping and the sum are checked.
        p = get_damping_param('pbe')
        for pair in ((18, 18), (18, 36)):
            q = np.asarray(data.sqrt_z_r4_over_r2)[list(pair)]  # indexed by atomic number
            rr = 3 * q[0] * q[1]
            f = p['a1'] * np.sqrt(rr) + p['a2']
            for distance in (3., 10., 50., 200.):  # 200 bohr is beyond dftd4's default 60 bohr cutoff
                with self.subTest(pair=pair, distance=distance):
                    coords = np.array([[0., 0, 0], [0, 0, distance]])
                    energy = d4_energy(pair, coords, 0, 'pbe')['energy_hartree']
                    c6 = DispersionModel(np.array(pair), coords).get_properties()['c6 coefficients'][0, 1]
                    expected = -c6 * (p['s6'] / (distance**6 + f**6) + p['s8'] * rr / (distance**8 + f**8))
                    self.assertLess(energy, 0)
                    self.assertAlmostEqual(energy / expected, 1., delta=1e-10)
        water = np.array([[0., 0, 0], [0, -1.43, 1.11], [0, 1.43, 1.11]])
        reference = d4_energy([8, 1, 1], water, 0, 'pbe')['energy_hartree']
        angle = .7
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        moved = water @ rotation.T + np.array([3., -2, 5])
        self.assertAlmostEqual(d4_energy([8, 1, 1], moved, 0, 'pbe')['energy_hartree'], reference, delta=1e-13)
        # EEQ charges, hence C6 and the energy, depend on the total charge.
        cation = d4_energy([8, 1, 1], water, 1, 'pbe')
        self.assertNotAlmostEqual(cation['energy_hartree'], reference, delta=1e-8)
        self.assertAlmostEqual(sum(cation['eeq_partial_charges']), 1., places=8)

    def test_scf_integration_both_engines_and_units(self):
        from engines import calculate
        from dftd4.pyscf import DFTD4Dispersion
        from pyscf import gto
        cfg = DFTConfig(atoms=parse_coordinates(WATER), basis='sto-3g', xc='pbe', grid_level=2,
                        cube_points=12, threads=1, conv_tol=1e-10)
        plain = calculate(cfg, self.root / 'plain', lambda _: None)
        d4 = calculate(replace(cfg, dispersion='d4'), self.root / 'd4', lambda _: None)
        native = calculate(replace(cfg, dispersion='d4', engine='native'), self.root / 'native', lambda _: None)
        # The dftd4 package's own PySCF adapter: independent of this project's unit/charge plumbing.
        mol = gto.M(atom=WATER, basis='sto-3g', verbose=0)
        adapter, _ = DFTD4Dispersion(mol, xc='pbe').kernel()
        for result in (d4, native):
            self.assertTrue(result['converged'])
            dispersion = result['dispersion']['energy_hartree']
            self.assertAlmostEqual(dispersion, adapter, delta=1e-12)
            self.assertLess(dispersion, 0)
            self.assertAlmostEqual(result['total_energy_hartree'], result['scf_energy_hartree'] + dispersion, delta=1e-12)
            self.assertAlmostEqual(result['energy_components_hartree']['dispersion'], dispersion, delta=1e-15)
            self.assertLess(abs(sum(result['energy_components_hartree'].values()) - result['total_energy_hartree']), 1e-9)
            self.assertAlmostEqual(result['scf_history'][-1]['energy_hartree'], result['scf_energy_hartree'], places=9)
        # Dispersion is post-SCF: the SCF energy is unchanged by switching it on.
        self.assertAlmostEqual(d4['scf_energy_hartree'], plain['total_energy_hartree'], delta=1e-9)
        self.assertIsNone(plain['dispersion'])
        self.assertLess(abs(native['total_energy_hartree'] - d4['total_energy_hartree']), 2e-7)
        bohr = replace(cfg, dispersion='d4', unit='bohr',
                       atoms=parse_coordinates('\n'.join(f'{a.element} ' + ' '.join(str(x / 0.52917721092) for x in a.position)
                                                         for a in cfg.atoms)))
        in_bohr = calculate(bohr, self.root / 'bohr', lambda _: None)
        self.assertAlmostEqual(in_bohr['dispersion']['energy_hartree'], d4['dispersion']['energy_hartree'], delta=1e-12)

    def test_missing_library_is_explicit_and_none_still_runs(self):
        from engines import calculate
        cfg = DFTConfig(atoms=parse_coordinates('H 0 0 0; H 0 0 .74'), basis='sto-3g', xc='pbe',
                        grid_level=1, cube_points=12, threads=1)
        blocked = {name: None for name in ('dftd4', 'dftd4.interface', 'dftd4.parameters')}
        with patch.dict(sys.modules, blocked):
            with self.assertRaisesRegex(RuntimeError, 'dftd4'):
                calculate(replace(cfg, dispersion='d4'), self.root / 'missing', lambda _: None)
            self.assertTrue(calculate(cfg, self.root / 'plain', lambda _: None)['converged'])


if __name__ == '__main__':
    unittest.main()
