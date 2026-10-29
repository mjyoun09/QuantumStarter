"""Per-element basis selection: input contract, checkpoint compatibility, real mixed-basis SCF."""
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from molecular import DFTConfig, describe_bases, parse_basis_overrides, parse_coordinates
from scf_control import signature

WATER = 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587'


class BasisOverrideInputTests(unittest.TestCase):
    def test_normalization_lookup_and_round_trip(self):
        cfg = DFTConfig(basis='sto-3g', basis_overrides={' o ': 'DEF2-SVP', 'H': 'sto-3g', 'Au': 'def2-svp'})
        # Same-as-default entries carry no physics and are dropped; keys are canonical symbols.
        self.assertEqual(cfg.basis_overrides, {'Au': 'def2-svp', 'O': 'def2-svp'})
        self.assertEqual(cfg.basis_by_element, {'H': 'sto-3g', 'O': 'def2-svp'})
        self.assertEqual(cfg.unused_basis_overrides, ['Au'])
        self.assertEqual(cfg.basis_label, 'sto-3g (O=def2-svp)')
        self.assertEqual(DFTConfig.from_dict(json.loads(json.dumps(cfg.to_dict()))), cfg)
        legacy = DFTConfig().to_dict()
        legacy.pop('basis_overrides')
        self.assertEqual(DFTConfig.from_dict(legacy), DFTConfig())
        self.assertEqual(describe_bases('sto-3g', {'H': 'sto-3g', 'O': 'def2-svp'}), 'sto-3g (O=def2-svp)')
        self.assertEqual(describe_bases('sto-3g', {}), 'sto-3g')

    def test_text_parser(self):
        self.assertEqual(parse_basis_overrides(' o=DEF2-TZVP; h = sto-3g\n'), {'O': 'def2-tzvp', 'H': 'sto-3g'})
        self.assertEqual(parse_basis_overrides(''), {})
        for text in ('O', 'O=a=b', 'O=def2-svp, o=sto-3g'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_basis_overrides(text)

    def test_rejections(self):
        cfg = DFTConfig()
        for overrides in ({'Xx': 'sto-3g'}, {'O': 'not-a-basis'}, {'O': 3}, {1: 'sto-3g'}):
            with self.subTest(overrides=overrides), self.assertRaises(ValueError):
                replace(cfg, basis_overrides=overrides)
        with self.assertRaises(ValueError):
            replace(cfg, basis_overrides=[('O', 'sto-3g')])
        # Og has no def2 data, even when the default basis would have been fine for it.
        with self.assertRaisesRegex(ValueError, 'Og'):
            DFTConfig(atoms=parse_coordinates('Og 0 0 0'), basis='dyall-v2z', basis_overrides={'Og': 'def2-svp'})
        # def2-ECP on Au cannot be combined with the X2C that Dyall on Pt switches on.
        with self.assertRaisesRegex(ValueError, 'Au'):
            DFTConfig(atoms=parse_coordinates('Au 0 0 0; Pt 0 0 2.5'), spin=1, basis='def2-svp',
                      basis_overrides={'Pt': 'dyall-v2z'})

    def test_ecp_relativity_and_grid_follow_each_element(self):
        heavy_ecp = DFTConfig(atoms=parse_coordinates('Au 0 0 0; H 0 0 1.52'), basis='sto-3g',
                              basis_overrides={'Au': 'def2-svp'})
        self.assertEqual(heavy_ecp.resolved_ecp, {'Au': 'def2-svp'})
        self.assertEqual((heavy_ecp.ecp_core_electrons, heavy_ecp.explicit_electrons), (60, 20))
        self.assertEqual(heavy_ecp.resolved_relativity, 'none')
        dyall = DFTConfig(atoms=parse_coordinates('Au 0 0 0; H 0 0 1.52'), basis='dyall-v2z',
                          basis_overrides={'H': 'def2-svp'})
        self.assertEqual(dyall.resolved_ecp, {})
        self.assertEqual((dyall.resolved_relativity, dyall.resolved_nuclear_model, dyall.resolved_grid_scheme),
                         ('x2c', 'gaussian', 'universal'))
        # A light element on def2 never receives an ECP; only the def2-ECP elements do.
        light = DFTConfig(basis='sto-3g', basis_overrides={'O': 'def2-svp'})
        self.assertEqual(light.resolved_ecp, {})
        with self.assertRaisesRegex(ValueError, 'def2-ECP'):
            replace(light, ecp='def2')

    def test_checkpoint_signature_compatibility(self):
        cfg = DFTConfig(basis='sto-3g')

        def legacy_signature(config):
            # The formula used before basis_overrides existed (scf_control at commit dd2871e),
            # applied to the settings that existed then: fields added later are removed first.
            data = config.to_dict()
            for key in ('engine', 'calculation', 'max_cycle', 'conv_tol', 'threads', 'memory_mb', 'cube_points',
                        'cube_margin_bohr', 'initial_guess', 'scf_solver', 'native_density_tol',
                        'native_residual_tol', 'native_damping', 'native_diis_space', 'opt_convergence', 'opt_max_steps', 'freq_step_bohr', 'thermo_temperature_k',
                        'basis_overrides', 'dispersion'):
                data.pop(key, None)
            return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()

        self.assertEqual(signature(cfg), legacy_signature(cfg))
        self.assertEqual(signature(replace(cfg, basis_overrides={'Au': 'def2-svp', 'H': 'sto-3g'})), signature(cfg))
        mixed = replace(cfg, basis_overrides={'O': 'def2-svp'})
        self.assertNotEqual(signature(mixed), signature(cfg))
        self.assertNotEqual(signature(mixed), signature(replace(cfg, basis_overrides={'H': 'def2-svp'})))


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class MixedBasisPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_matches_direct_pyscf_and_native_engine(self):
        from pyscf import dft, gto, lib
        from pyscf.dft import rks
        from engines import calculate
        from engines.pyscf_engine import configure_grids
        lib.num_threads(1)
        cfg = DFTConfig(atoms=parse_coordinates(WATER), basis='sto-3g', basis_overrides={'O': 'def2-svp'},
                        xc='lda,vwn', grid_level=2, cube_points=12, threads=1, conv_tol=1e-10)
        external = calculate(cfg, self.root / 'pyscf', lambda _: None)
        with patch.object(rks.RKS, 'kernel', side_effect=AssertionError('External SCF called')):
            native = calculate(replace(cfg, engine='native'), self.root / 'native', lambda _: None)
        uniform = calculate(replace(cfg, basis_overrides={}), self.root / 'uniform', lambda _: None)

        # Independent construction with PySCF's own dict-basis interface.
        mol = gto.M(atom=WATER, basis={'O': 'def2-svp', 'H': 'sto-3g'}, verbose=0)
        mf = dft.RKS(mol)
        mf.xc, mf.small_rho_cutoff, mf.conv_tol = 'lda,vwn', 0, 1e-10
        configure_grids(mf.grids, cfg)
        direct = mf.kernel()
        self.assertTrue(mf.converged)
        # def2-SVP O: 3s2p1d = 14 spherical AO; STO-3G H: 1 each.
        self.assertEqual((mol.nao_nr(), external['basis_functions'], native['basis_functions']), (16, 16, 16))
        self.assertEqual(uniform['basis_functions'], 7)
        for result in (external, native):
            self.assertTrue(result['converged'])
            self.assertEqual(result['basis_by_element'], {'H': 'sto-3g', 'O': 'def2-svp'})
            self.assertAlmostEqual(result['electron_count_from_density'], 10, places=8)
            self.assertTrue(any('원소별 기저 혼합' in w for w in result['warnings']))
        self.assertLess(abs(external['total_energy_hartree'] - direct), 1e-8)
        self.assertLess(abs(native['total_energy_hartree'] - external['total_energy_hartree']), 2e-7)
        # More variational freedom on O lowers the energy for the same functional and grid.
        self.assertLess(external['total_energy_hartree'], uniform['total_energy_hartree'])

    def test_heavy_element_ecp_on_one_element_only(self):
        from engines import calculate
        cfg = DFTConfig(atoms=parse_coordinates('Au 0 0 0; H 0 0 1.52'), basis='sto-3g',
                        basis_overrides={'Au': 'def2-svp'}, xc='lda,vwn', grid_level=2, cube_points=12, threads=1)
        result = calculate(cfg, self.root / 'auh', lambda _: None)
        self.assertTrue(result['converged'])
        self.assertEqual(result['ecp']['resolved'], {'Au': 'def2-svp'})
        self.assertEqual(result['ecp']['core_electrons_by_atom'], [60, 0])
        self.assertEqual(result['explicit_electrons'], 20)
        self.assertAlmostEqual(result['electron_count_from_density'], 20, places=8)


if __name__ == '__main__':
    unittest.main()
