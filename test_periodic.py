"""Boundary geometry, k-point spin conventions, real PBC SCF and provenance."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from ase import Atoms
from ase.build import bulk, fcc111
from ase.io import write
from ase.calculators.calculator import Calculator, all_changes

from periodic import PeriodicAtom, PeriodicConfig, load_calculation_config, read_periodic_structure
from molecular import save_dft_config
from surface_workflow import (generate_candidates, read_json, relative_comparison, screen_candidates,
                              select_candidates, structure_hash)


class PeriodicInputTests(unittest.TestCase):
    def setUp(self):
        self.config = PeriodicConfig.from_ase(bulk('Si', 'diamond', a=5.43))

    def test_element_input_and_periodic_data_are_distinct(self):
        atom = PeriodicAtom('Og', (0, 0, 0))
        self.assertEqual(atom.element, 'Og')
        with self.assertRaisesRegex(ValueError, '기저|데이터'):
            replace(self.config, atoms=(atom,))

    def test_lattice_pbc_and_duplicate_images(self):
        for changes in (dict(lattice=np.zeros((3, 3))), dict(lattice=np.eye(3) * -3),
                        dict(lattice=np.eye(3) * np.nan), dict(kmesh=(1, 0, 1)),
                        dict(charge=1), dict(smearing_ev=-.1)):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                replace(self.config, **changes)
        with self.assertRaisesRegex(ValueError, '중복'):
            PeriodicConfig(atoms=(PeriodicAtom('Pt', (0, 0, 0)), PeriodicAtom('Pt', (3, 0, 0))), lattice=np.eye(3) * 3)
        slab = PeriodicConfig.from_ase(fcc111('Pt', (1, 1, 2), a=3.92, vacuum=8))
        self.assertEqual(slab.pbc, (True, True, False))
        with self.assertRaisesRegex(ValueError, 'z 방향'):
            replace(slab, kmesh=(2, 2, 2))
        with self.assertRaises(ValueError):
            replace(slab, lattice=((3, 0, 1), (0, 3, 0), (0, 0, 20)))

    def test_structure_and_config_roundtrip(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            save_dft_config(self.config, p / 'settings.json')
            self.assertEqual(load_calculation_config(p / 'settings.json'), self.config)
            atoms = bulk('Au', 'fcc', a=4.08)
            write(str(p / 'POSCAR'), atoms, format='vasp')
            loaded = read_periodic_structure(p / 'POSCAR', spin=1)
            self.assertEqual(loaded.atoms[0].element, 'Au')
            np.testing.assert_allclose(loaded.lattice, atoms.cell.array)

    def test_site_generator_is_deterministic_and_same_composition(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            first = read_json(generate_candidates(p / 'one', seed=0))
            second = read_json(generate_candidates(p / 'two', seed=0))
            self.assertEqual(first, second)
            self.assertEqual(len(first['candidates']), 36)
            self.assertEqual(first['substrate_composition'], dict(Ag=8, Au=7, Cu=7, Pd=7, Pt=7))

    def test_selection_excludes_failure_and_explores(self):
        rows = [dict(id=str(i), energy_ev=i, status='converged') for i in range(8)]
        rows += [dict(id='failed', energy_ev=-100, status='not_converged')]
        selected = select_candidates(rows, 2, 1, 42)
        self.assertEqual(selected[:2], ['0', '1'])
        self.assertEqual(len(selected), 3)
        self.assertNotIn('failed', selected)
        rows[1]['equivalent_minimum_of'] = '0'
        self.assertEqual(select_candidates(rows, 2, 0), ['0', '2'])

    def test_comparison_uses_energy_differences_and_subset(self):
        rows = [dict(id='a', ml_energy_ev=-1000, dft_energy_ev=-10, status='converged'),
                dict(id='b', ml_energy_ev=-999.7, dft_energy_ev=-9.8, status='converged'),
                dict(id='bad', ml_energy_ev=-1001, dft_energy_ev=-20, status='not_converged')]
        comparison = relative_comparison(rows)
        self.assertEqual(comparison['count'], 2)
        self.assertAlmostEqual(comparison['relative_energy_mae_ev'], .05)
        self.assertEqual(comparison['pair_order_agreement'], 1.)

    def test_unconverged_relaxations_and_changed_inputs_not_ranked(self):
        class ConstantForce(Calculator):
            implemented_properties = ['energy', 'forces']
            def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
                super().calculate(atoms, properties, system_changes)
                self.results = dict(energy=0., forces=np.ones((len(atoms), 3)))
        with tempfile.TemporaryDirectory() as temp:
            manifest = generate_candidates(Path(temp) / 'candidates', material='Pt', size=(1, 1, 2))
            result = read_json(screen_candidates(manifest, steps=0, calculator=(ConstantForce(), {'model': 'test-only'})))
            self.assertEqual(result['selected_ids'], [])
            self.assertTrue(all(r['status'] == 'not_converged' for r in result['rows']))
            item = read_json(manifest)['candidates'][0]
            path = manifest.parent / item['file']
            from ase.io import read
            changed = read(str(path), format='extxyz')
            changed.positions[-1, 2] += .01
            write(str(path), changed, format='extxyz')
            result = read_json(screen_candidates(manifest, steps=0, calculator=(ConstantForce(), {'model': 'test-only'})))
            self.assertEqual(result['rows'][0]['status'], 'error')
            self.assertIn('변경', result['rows'][0]['error'])


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires Linux PySCF')
class PeriodicEngineTests(unittest.TestCase):
    def test_target_metals_have_matching_valence_data(self):
        from engines.periodic import build_cell
        for symbol, expected in [('Cu', 11), ('Ag', 11), ('Au', 11), ('Pd', 18), ('Pt', 18)]:
            cfg = PeriodicConfig.from_ase(bulk(symbol, 'fcc', a=4), spin=expected % 2)
            cell, valence = build_cell(cfg)
            self.assertEqual(valence[symbol], expected)
            self.assertEqual(cell.nelectron, expected)
        with self.assertRaisesRegex(ValueError, '홀짝'):
            build_cell(PeriodicConfig.from_ase(bulk('Au', 'fcc', a=4)))

    def test_kpoint_spin_is_per_cell_and_nonconvergence_is_saved(self):
        from engines.periodic import calculate_periodic
        atoms = Atoms('H', positions=[[1, 1, 1]], cell=np.eye(3) * 4, pbc=True)
        cfg = PeriodicConfig.from_ase(atoms, spin=1, kmesh=(2, 1, 1), basis='gth-dzvp',
                                      grid_level=0, precision=1e-5, conv_tol=1e-6, max_cycle=50, threads=2)
        with tempfile.TemporaryDirectory() as temp:
            result = calculate_periodic(cfg, Path(temp) / 'spin', progress=lambda _: None)
            self.assertTrue(result['converged'])
            self.assertAlmostEqual(result['electron_trace'], 1, places=6)
            self.assertAlmostEqual(result['spin_trace_per_cell'], 1, places=6)
            self.assertEqual(len(result['kpoints_scaled']), 2)
            result = calculate_periodic(replace(cfg, max_cycle=1), Path(temp) / 'short', progress=lambda _: None)
            self.assertFalse(result['converged'])
            self.assertEqual(read_json(Path(temp) / 'short/results.json')['status'], 'not_converged')


if __name__ == '__main__':
    unittest.main()
