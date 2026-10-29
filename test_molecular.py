"""Input boundaries, physical invariants and real PySCF integration checks."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from dft_workflow import CalculationError, list_dft_runs, run_calculation
from molecular import (Atom, ATOMIC_NUMBERS, BASES, DEF2_CORE_ELECTRONS, DFTConfig, PERIOD4_ELEMENTS, SUPPORTED_ELEMENTS,
                       available_bases, load_dft_config, parse_coordinates, read_xyz,
                       save_dft_config)
from molecular_plotting import BOHR_TO_ANGSTROM, create_molecular_figure, read_cube


class MolecularInputTests(unittest.TestCase):
    def test_electron_charge_spin_and_method(self):
        water = DFTConfig()
        self.assertEqual(water.electrons, 10)
        self.assertEqual(water.resolved_method, 'rks')
        radical = replace(water, charge=1, spin=1)
        self.assertEqual(radical.electrons, 9)
        self.assertEqual(radical.resolved_method, 'uks')
        for changes in ({'spin': 1}, {'spin': -1}, {'spin': 12}, {'charge': 11},
                        {'charge': True}, {'spin': 2, 'method': 'rks'}, {'max_cycle': 0},
                        {'conv_tol': float('nan')}, {'unit': 'nm'}, {'basis': 'fake'},
                        {'xc': 'HF'}, {'system': 'periodic'}, {'scf_solver': 'fake'}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                replace(water, **changes)

    def test_atom_and_coordinate_validation(self):
        self.assertEqual(parse_coordinates('he 0 0 0')[0].element, 'He')
        for text in ('', 'H 0 0', 'H zero 0 0', 'H nan 0 0', 'Qq 0 0 0', 'Uue 0 0 0'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_coordinates(text)
        with self.assertRaisesRegex(ValueError, '중복'):
            DFTConfig(atoms=parse_coordinates('H 0 0 0; H 0 0 0'))

    def test_period4_boundary_and_charge_spin(self):
        self.assertEqual(PERIOD4_ELEMENTS, ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
                                           'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'))
        for z, symbol in enumerate(PERIOD4_ELEMENTS, 19):
            with self.subTest(element=symbol):
                cfg = DFTConfig(atoms=parse_coordinates(f'{symbol} 0 0 0'), spin=z % 2)
                self.assertEqual(cfg.electrons, z)
        fe = DFTConfig(atoms=parse_coordinates('Fe 0 0 0'), charge=3, spin=5)
        co = DFTConfig(atoms=parse_coordinates('Co 0 0 0'), charge=2, spin=3)
        self.assertEqual((fe.electrons, co.electrons), (23, 25))
        self.assertEqual(co.resolved_method, 'uks')
        for cfg in (fe, co):
            self.assertEqual(DFTConfig.from_dict(cfg.to_dict()), cfg)
            with self.assertRaises(ValueError):
                replace(cfg, spin=0)

    def test_all_118_elements_and_relativistic_configuration(self):
        self.assertEqual(len(SUPPORTED_ELEMENTS), 118)
        self.assertEqual(SUPPORTED_ELEMENTS[-1], 'Og')
        for z, symbol in enumerate(SUPPORTED_ELEMENTS, 1):
            with self.subTest(element=symbol):
                cfg = DFTConfig(atoms=parse_coordinates(f'{symbol} 0 0 0'), spin=z % 2, basis='dyall-v2z')
                self.assertEqual(cfg.explicit_electrons, z)
                self.assertEqual(cfg.resolved_relativity, 'x2c')
                self.assertEqual(cfg.resolved_nuclear_model, 'gaussian')
                self.assertEqual(cfg.resolved_grid_scheme, 'universal')
                self.assertEqual(DFTConfig.from_dict(cfg.to_dict()), cfg)
        self.assertEqual(available_bases(cfg.atoms), ('dyall-v2z', 'dyall-v3z'))
        for changes in ({'basis': 'def2-svp'}, {'grid_scheme': 'treutler'}, {'relativity': 'soc'},
                        {'nuclear_model': 'fake'}, {'grid_scheme': 'fake'}):
            with self.assertRaises(ValueError):
                replace(cfg, **changes)
        for symbol, core in (('Cs', 46), ('La', 46), ('Hf', 60), ('Au', 60), ('Rn', 60)):
            cfg = DFTConfig(atoms=parse_coordinates(f'{symbol} 0 0 0'), spin=ATOMIC_NUMBERS[symbol] % 2)
            self.assertEqual(cfg.ecp_core_electrons, core)
            for changes in ({'ecp': 'none'}, {'relativity': 'x2c'}, {'nuclear_model': 'gaussian'}):
                with self.assertRaises(ValueError):
                    replace(cfg, **changes)
        for symbol in ('Ce', 'Lu', 'U', 'Og'):
            with self.assertRaisesRegex(ValueError, '기저 데이터'):
                DFTConfig(atoms=parse_coordinates(f'{symbol} 0 0 0'))

    def test_ecp_pairing_explicit_electrons_and_roundtrip(self):
        for symbol, z in (('Rb', 37), ('Sr', 38), ('Y', 39), ('Zr', 40)):
            cfg = DFTConfig(atoms=parse_coordinates(f'{symbol} 0 0 0'), spin=z % 2)
            self.assertEqual((cfg.electrons, cfg.ecp_core_electrons, cfg.explicit_electrons), (z, 28, z - 28))
            self.assertEqual(cfg.resolved_ecp, {symbol: 'def2-svp'})
            self.assertEqual(DFTConfig.from_dict(cfg.to_dict()), cfg)
            self.assertEqual(replace(cfg, ecp='def2').resolved_ecp, cfg.resolved_ecp)
            with self.assertRaisesRegex(ValueError, 'ECP'):
                replace(cfg, ecp='none')
            full = replace(cfg, basis='sto-3g', ecp='none')
            self.assertEqual(full.explicit_electrons, z)
            self.assertEqual(full.resolved_ecp, {})
        zr = DFTConfig(atoms=parse_coordinates('Zr 0 0 0'), charge=4)
        self.assertEqual((zr.electrons, zr.explicit_electrons), (36, 8))
        for changes in ({'spin': 10}, {'charge': 12}, {'ecp': 'lanl2dz'}, {'basis': 'sto-3g', 'ecp': 'def2'}):
            with self.assertRaises(ValueError):
                replace(zr, **changes)
        with self.assertRaises(ValueError):
            DFTConfig(ecp='def2')
        mixed = DFTConfig(atoms=parse_coordinates('Zr 0 0 0; O 0 0 2; Zr 0 0 4'))
        self.assertEqual((mixed.electrons, mixed.ecp_core_electrons, mixed.explicit_electrons), (88, 56, 32))
        self.assertEqual(available_bases(mixed.atoms), ('sto-3g', 'def2-svp', 'def2-tzvp', 'dyall-v2z', 'dyall-v3z'))
        with tempfile.TemporaryDirectory() as temp:
            file = Path(temp) / 'ecp.json'
            save_dft_config(mixed, file)
            self.assertEqual(load_dft_config(file), mixed)

    def test_basis_compatibility_for_mixed_molecules(self):
        kcl = parse_coordinates('K 0 0 0; Cl 0 0 2.67')
        self.assertNotIn('cc-pvdz', available_bases(kcl))
        with self.assertRaisesRegex(ValueError, 'cc-pvdz.*K'):
            DFTConfig(atoms=kcl, basis='cc-pvdz')
        for symbol in ('Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co'):
            atoms = parse_coordinates(f'{symbol} 0 0 0')
            self.assertNotIn('6-31+g*', available_bases(atoms))
            self.assertIn('def2-svp', available_bases(atoms))
            self.assertIn('def2-tzvp', available_bases(atoms))
            with self.assertRaisesRegex(ValueError, f'6-31\\+g\\*.*{symbol}'):
                DFTConfig(atoms=atoms, basis='6-31+g*')
        self.assertEqual(available_bases(DFTConfig().atoms), BASES)

    def test_xyz_frame_and_config_roundtrip(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / '구조.xyz'
            path.write_text('2\n\nH 0 0 0\nH 0 0 .74\n', encoding='utf-8-sig')
            cfg = DFTConfig(atoms=read_xyz(path))
            self.assertEqual(cfg.electrons, 2)
            settings = Path(temp) / 'settings.json'
            save_dft_config(cfg, settings)
            self.assertEqual(load_dft_config(settings), cfg)
            for text in ('3\ncomment\nH 0 0 0\n', '1\ncomment\nH 0 0 0\n1\nnext\nH 0 0 1',
                         '2\ncomment\nH 0 0 0\n\n'):
                path.write_text(text, encoding='utf-8')
                with self.assertRaises(ValueError):
                    read_xyz(path)
            with self.assertRaises(ValueError):
                DFTConfig.from_dict({'nuclear_charge': 1})
            with self.assertRaises(ValueError):
                DFTConfig.from_dict({'atoms': [{'element': 'H', 'position': [0, 0]}]})

    def test_failed_backend_keeps_settings_and_failed_history(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch('dft_workflow.worker_command', side_effect=RuntimeError('엔진 설치 필요')):
                with self.assertRaises(CalculationError) as error:
                    run_calculation(DFTConfig(), temp)
            folder = error.exception.directory
            self.assertTrue((folder / 'settings.json').is_file())
            self.assertFalse((folder / 'results.json').exists())
            self.assertEqual(list_dft_runs(temp)[0]['summary']['status'], 'error')


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'PySCF integration runs in .venv-dft (WSL)')
class MolecularPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def calculate(self, config, name):
        from engines import calculate
        return calculate(config, self.root / name, lambda _: None)

    def check_electronic_invariants(self, summary, name):
        self.assertTrue(summary['converged'])
        data = np.load(self.root / name / 'electronic_structure.npz', allow_pickle=False)
        dm, coeff, overlap = data['density_matrix'], data['mo_coeff'], data['overlap']
        occ = data['mo_occupation']
        self.assertAlmostEqual(float(occ.sum()), summary['explicit_electrons'], places=8)
        self.assertAlmostEqual(summary['electron_count_from_density'], summary['explicit_electrons'], places=8)
        self.assertEqual(summary['explicit_electrons'] + summary['ecp_core_electrons'], summary['electrons'])
        self.assertAlmostEqual(summary['spin_from_density'], summary['config']['spin'], places=8)
        coeffs = [coeff] if coeff.ndim == 2 else list(coeff)
        dms = [dm] if dm.ndim == 2 else list(dm)
        occs = [occ] if occ.ndim == 1 else list(occ)
        for c, d, o in zip(coeffs, dms, occs):
            np.testing.assert_allclose(c.T @ overlap @ c, np.eye(c.shape[1]), atol=1e-8)
            np.testing.assert_allclose((c * o) @ c.T, d, atol=1e-8)
        if dm.ndim == 3:
            self.assertAlmostEqual(float(np.einsum('ij,ji->', dm[0], overlap)), summary['alpha_electrons'], places=8)
            self.assertAlmostEqual(float(np.einsum('ij,ji->', dm[1], overlap)), summary['beta_electrons'], places=8)
        density, _, _ = read_cube(self.root / name / 'density.cube')
        self.assertGreaterEqual(float(density.min()), -1e-10)
        self.assertEqual(density.shape, (16, 16, 16))
        self.assertEqual(load_dft_config(self.root / name / 'settings.json'), DFTConfig.from_dict(summary['config']))
        self.assertTrue(all((self.root / name / file).is_file() for file in summary['cube_files']))
        numbers = [ATOMIC_NUMBERS[a['element']] for a in summary['config']['atoms']]
        np.testing.assert_array_equal(data['atom_numbers'], numbers)
        np.testing.assert_array_equal(data['effective_nuclear_charges'] + data['ecp_core_electrons_by_atom'], numbers)
        with np.load(self.root / name / 'density_grid.npz', allow_pickle=False) as grid:
            np.testing.assert_array_equal(grid['atom_numbers'], numbers)
        for file in summary['cube_files']:
            lines = (self.root / name / file).read_text().splitlines()
            self.assertEqual([int(lines[6 + i].split()[0]) for i in range(len(numbers))], numbers)
            if summary['ecp_core_electrons']:
                self.assertIn('core electrons omitted', lines[0])

    def test_h2_rks_energy_density_and_basis_functional_replacement(self):
        cfg = DFTConfig(atoms=parse_coordinates('H 0 0 0; H 0 0 .74'), basis='sto-3g',
                        grid_level=3, cube_points=16, threads=1)
        summary = self.calculate(cfg, 'h2')
        # PBE/STO-3G H2 at 0.74 Å; loose across PySCF/Libxc versions.
        self.assertAlmostEqual(summary['total_energy_hartree'], -1.1520728, delta=2e-5)
        self.check_electronic_invariants(summary, 'h2')
        hybrid = self.calculate(replace(cfg, basis='def2-svp', xc='b3lyp'), 'hybrid')
        self.check_electronic_invariants(hybrid, 'hybrid')
        self.assertGreater(hybrid['basis_functions'], summary['basis_functions'])
        self.assertGreater(abs(hybrid['total_energy_hartree'] - summary['total_energy_hartree']), .01)

    def test_neutral_helium_open_shell_hydrogen_and_ion(self):
        for name, text, charge, spin, electrons in (
            ('he', 'He 0 0 0', 0, 0, 2), ('h', 'H 0 0 0', 0, 1, 1),
            ('h2plus', 'H 0 0 0; H 0 0 1.1', 1, 1, 1)):
            cfg = DFTConfig(atoms=parse_coordinates(text), charge=charge, spin=spin,
                            basis='sto-3g', cube_points=16, threads=1)
            summary = self.calculate(cfg, name)
            self.check_electronic_invariants(summary, name)
            self.assertEqual(summary['electrons'], electrons)
            if spin:
                self.assertEqual(summary['method'], 'UKS')
                self.assertTrue((self.root / name / 'spin_density.cube').exists())
                self.assertAlmostEqual(summary['spin_square'], .75, places=7)

    def test_water_units_and_translation_invariance(self):
        cfg = DFTConfig(basis='sto-3g', cube_points=16, threads=1)
        a = self.calculate(cfg, 'water')
        self.check_electronic_invariants(a, 'water')
        shifted = tuple(Atom(atom.element, tuple((x + offset) / BOHR_TO_ANGSTROM
                        for x, offset in zip(atom.position, (1.2, -.4, 2)))) for atom in cfg.atoms)
        b = self.calculate(replace(cfg, atoms=shifted, unit='bohr'), 'shifted')
        self.assertAlmostEqual(a['total_energy_hartree'], b['total_energy_hartree'], delta=1e-7)
        for field in a['cube_files']:
            for plane in ('xy', 'xz', 'yz'):
                fig = create_molecular_figure(a, self.root / 'water', field, plane)
                self.assertEqual(len(fig.axes), 3)

    def test_nonconvergence_never_marked_success_and_cli_exit_code(self):
        cfg = DFTConfig(basis='sto-3g', max_cycle=1, conv_tol=1e-12, cube_points=16, threads=1)
        summary, directory = run_calculation(cfg, self.root, backend='local')
        self.assertFalse(summary['converged'])
        self.assertEqual(summary['status'], 'not_converged')
        self.assertEqual(list_dft_runs(self.root)[0]['summary']['status'], 'not_converged')
        self.assertTrue((directory / 'density.cube').is_file())
        self.assertTrue(any('미수렴' in w for w in summary['warnings']))
        self.assertFalse((directory / 'error.json').exists())

    def test_supported_basis_element_contract(self):
        from pyscf import gto
        import warnings
        supported = 0
        # Check the app's capability table against actual installed basis data.
        for symbol in SUPPORTED_ELEMENTS:
            atom = Atom(symbol, (0, 0, 0))
            for basis in BASES:
                with self.subTest(element=symbol, basis=basis), warnings.catch_warnings():
                    warnings.simplefilter('ignore', UserWarning)
                    try:
                        loaded = bool(gto.basis.load(basis, symbol))
                    except gto.basis.BasisNotFoundError:
                        loaded = False
                    self.assertEqual(loaded, basis in available_bases((atom,)))
                    supported += int(loaded)
        self.assertEqual(supported, 554)

    def test_all_element_molecule_and_grid_construction(self):
        from pyscf import dft, gto, lib
        from engines.pyscf_engine import configure_grids
        lib.num_threads(1)
        for z, symbol in enumerate(SUPPORTED_ELEMENTS, 1):
            for basis in ('dyall-v2z', 'dyall-v3z'):
                with self.subTest(element=symbol, basis=basis):
                    cfg = DFTConfig(atoms=(Atom(symbol, (0, 0, 0)),), spin=z % 2, basis=basis, grid_level=0)
                    mol = gto.M(atom=[(symbol, (0, 0, 0))], spin=cfg.spin, basis=basis, nucmod='gaussian', verbose=0)
                    self.assertEqual(mol.nelectron, z)
                    self.assertLessEqual(mol.nao_nr(), 800)
                    self.assertGreater(mol._env[mol._atm[0, gto.PTR_ZETA]], 0)
                    self.assertTrue(np.isfinite(mol.intor_symmetric('int1e_ovlp')).all())
                    if basis == 'dyall-v2z':
                        grids = dft.gen_grid.Grids(mol)
                        configure_grids(grids, cfg)
                        grids.build()
                        self.assertGreater(len(grids.weights), 0)
                        self.assertTrue(np.isfinite(grids.weights).all())
                        self.assertTrue(np.isfinite(grids.coords).all())

    def test_x2c_closed_open_shell_and_nonconvergence(self):
        for symbol, spin, solver in (('He', 0, 'diis'), ('H', 1, 'newton')):
            cfg = DFTConfig(atoms=(Atom(symbol, (0, 0, 0)),), spin=spin, basis='dyall-v2z',
                            scf_solver=solver, grid_level=2, cube_points=16, threads=1)
            result = self.calculate(cfg, symbol)
            self.check_electronic_invariants(result, symbol)
            self.assertEqual(result['relativity'], 'x2c')
            self.assertEqual(result['nuclear_model'], 'gaussian')
            self.assertEqual(result['integration_grid']['scheme'], 'universal')
            self.assertLess(abs(result['integration_grid']['electron_count_error']), .01)
            self.assertEqual(result['density']['representation'], 'x2c_picture_density')
            self.assertIn('X2C', (self.root / symbol / 'density.cube').read_text().splitlines()[0])
            self.assertIn('X2C', create_molecular_figure(result, self.root / symbol).axes[1].get_title())
        cfg = DFTConfig(basis='dyall-v2z', max_cycle=1, conv_tol=1e-12, cube_points=16, grid_level=1, threads=1)
        result = self.calculate(cfg, 'incomplete_x2c')
        self.assertFalse(result['converged'])
        self.assertEqual(result['status'], 'not_converged')
        self.assertAlmostEqual(result['electron_count_from_density'], 10, places=8)

    def test_def2_ecp_data_contract(self):
        from pyscf import gto
        for basis in ('def2-svp', 'def2-tzvp'):
            for symbol, core in DEF2_CORE_ELECTRONS.items():
                with self.subTest(basis=basis, element=symbol):
                    cfg = DFTConfig(atoms=parse_coordinates(f'{symbol} 0 0 0'),
                                    spin=ATOMIC_NUMBERS[symbol] % 2, basis=basis)
                    self.assertEqual(gto.basis.load_ecp(basis, symbol)[0], core)
                    mol = gto.M(atom=[(symbol, (0, 0, 0))], basis=basis,
                                ecp=cfg.resolved_ecp, spin=cfg.spin, verbose=0)
                    self.assertEqual(mol.nelectron, cfg.explicit_electrons)
                    self.assertEqual(mol.atom_nelec_core(0), core)

    def test_elements_28_to_40_and_ecp_physics(self):
        # Chosen charge/spin states exercise each added element, including open-shell ECP.
        # These are numerical invariants, not a claim of accurate ground-state energies.
        cases = [('Ni', 2, 2), ('Cu', 1, 0), ('Zn', 2, 0), ('Ga', 3, 0), ('Ge', 4, 0),
                 ('As', 3, 0), ('Se', 4, 0), ('Br', -1, 0), ('Kr', 0, 0),
                 ('Rb', 0, 1), ('Sr', 0, 0), ('Y', 3, 0), ('Zr', 4, 0)]
        for symbol, charge, spin in cases:
            with self.subTest(element=symbol, charge=charge):
                cfg = DFTConfig(atoms=parse_coordinates(f'{symbol} 0 0 0'), charge=charge, spin=spin,
                                basis='def2-svp', cube_points=16, grid_level=2, threads=1,
                                scf_solver='newton' if symbol == 'Ni' else 'diis')
                summary = self.calculate(cfg, symbol)
                self.check_electronic_invariants(summary, symbol)
                self.assertEqual(summary['electrons'], ATOMIC_NUMBERS[symbol] - charge)
                self.assertEqual(summary['ecp_core_electrons'], 28 if symbol in ('Rb', 'Sr', 'Y', 'Zr') else 0)
                if cfg.resolved_ecp:
                    self.assertFalse(summary['density']['includes_core_electrons'])
                    self.assertEqual(summary['density']['representation'], 'ecp_pseudo_density')
                    self.assertTrue(any('def2-ECP' in text for text in summary['warnings']))
                if symbol == 'Zr':
                    figure = create_molecular_figure(summary, self.root / symbol)
                    self.assertIn('ECP', figure.axes[1].get_title())
                    self.assertTrue(any(t.get_text().strip() == 'Zr' for t in figure.axes[0].texts))

    def test_mixed_ecp_molecule_and_newton_nonconvergence(self):
        cfg = DFTConfig(atoms=parse_coordinates('Rb 0 0 0; Cl 0 0 2.8'), cube_points=16,
                        grid_level=2, threads=1)
        summary = self.calculate(cfg, 'rbcl')
        self.check_electronic_invariants(summary, 'rbcl')
        self.assertEqual((summary['electrons'], summary['explicit_electrons']), (54, 26))
        self.assertEqual(summary['ecp']['core_electrons_by_atom'], [28, 0])
        self.assertEqual(summary['ecp']['effective_nuclear_charges'], [9, 17])
        incomplete = self.calculate(replace(cfg, scf_solver='newton', max_cycle=1, conv_tol=1e-12), 'rbcl_incomplete')
        self.assertFalse(incomplete['converged'])
        self.assertEqual(incomplete['status'], 'not_converged')
        self.assertAlmostEqual(incomplete['electron_count_from_density'], 26, places=8)

    def test_period4_atoms_and_transition_metal_ions(self):
        for name, symbol, charge, spin, electrons in (
            ('potassium', 'K', 0, 1, 19), ('calcium', 'Ca', 0, 0, 20),
            ('ferric', 'Fe', 3, 5, 23), ('cobalt', 'Co', 2, 3, 25)):
            with self.subTest(element=symbol, charge=charge):
                cfg = DFTConfig(atoms=parse_coordinates(f'{symbol} 0 0 0'), charge=charge, spin=spin,
                                basis='def2-svp', cube_points=16, grid_level=2, threads=1,
                                scf_solver='newton' if symbol == 'Co' else 'diis')
                summary = self.calculate(cfg, name)
                self.check_electronic_invariants(summary, name)
                self.assertEqual(summary['electrons'], electrons)
                self.assertEqual(summary['target_spin_square'], spin / 2 * (spin / 2 + 1))
                if symbol in ('Fe', 'Co'):
                    self.assertTrue(any('전이금속' in text for text in summary['warnings']))
                    self.assertTrue((self.root / name / 'spin_density.cube').exists())
                if symbol == 'Co':
                    self.assertEqual(summary['scf_solver'], 'newton')
                    self.assertEqual(summary['scf_iteration_kind'], 'macroiteration')
                    self.assertEqual([r['cycle'] for r in summary['scf_history']],
                                     list(range(1, summary['scf_cycles'] + 1)))
                    self.assertTrue(all(r['density_change_norm'] is None for r in summary['scf_history']))

    def test_newton_nonconvergence_preserved(self):
        cfg = DFTConfig(atoms=parse_coordinates('Co 0 0 0'), charge=2, spin=3, basis='def2-svp',
                        scf_solver='newton', max_cycle=1, conv_tol=1e-12,
                        cube_points=16, grid_level=2, threads=1)
        summary = self.calculate(cfg, 'newton_incomplete')
        self.assertFalse(summary['converged'])
        self.assertEqual(summary['status'], 'not_converged')
        self.assertEqual(summary['scf_cycles'], 1)
        self.assertEqual(summary['scf_solver'], 'newton')
        self.assertAlmostEqual(summary['spin_from_density'], 3, places=8)

    def test_kcl_molecular_calculation(self):
        cfg = DFTConfig(atoms=parse_coordinates('K 0 0 0; Cl 0 0 2.67'), basis='def2-svp',
                        cube_points=16, grid_level=2, threads=1)
        summary = self.calculate(cfg, 'kcl')
        self.check_electronic_invariants(summary, 'kcl')
        self.assertEqual(summary['electrons'], 36)
        self.assertEqual(summary['method'], 'RKS')
        self.assertGreater(summary['nuclear_repulsion_hartree'], 0)


if __name__ == '__main__':
    unittest.main()
