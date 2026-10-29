"""Harmonic frequencies: projection/units vs PySCF, analytic Hessian, H2 curvature, saddle point."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from molecular import DFTConfig, parse_coordinates
from scf_control import signature

WATER = 'O 0 0 0; H 0 -.80 .55; H 0 .80 .55'
READY = bool(importlib.util.find_spec('pyscf') and importlib.util.find_spec('geometric'))


class FrequencyInputTests(unittest.TestCase):
    def test_fields(self):
        cfg = DFTConfig(calculation='frequencies', freq_step_bohr=0.01, thermo_temperature_k=500)
        self.assertEqual(DFTConfig.from_dict(json.loads(json.dumps(cfg.to_dict()))), cfg)
        self.assertEqual(signature(cfg), signature(DFTConfig()))
        for kwargs in ({'freq_step_bohr': 0}, {'freq_step_bohr': 0.1}, {'thermo_temperature_k': 0}):
            with self.subTest(**kwargs), self.assertRaises(ValueError):
                replace(cfg, **kwargs)

    def test_thermo_formulas(self):
        from engines.frequencies import CM1_TO_HARTREE, vibrational_thermo
        result = vibrational_thermo([-50., 1000., 2000.], 298.15)
        self.assertAlmostEqual(result['zpe_hartree'], 1500 * CM1_TO_HARTREE, places=15)
        self.assertGreater(result['u_vib_hartree'], result['zpe_hartree'])
        self.assertGreater(result['s_vib_hartree_per_k'], 0)
        cold = vibrational_thermo([1000.], 1.0)
        self.assertAlmostEqual(cold['u_vib_hartree'], cold['zpe_hartree'], places=15)
        self.assertTrue(all(np.isfinite(v) for v in cold.values() if isinstance(v, float)))
        self.assertEqual((cold['s_vib_hartree_per_k'], cold['cv_vib_hartree_per_k']), (0.0, 0.0))
        # High-temperature limit: each mode's heat capacity approaches k_B.
        from scipy.constants import Boltzmann, physical_constants
        hot = vibrational_thermo([10.], 5000.)
        self.assertAlmostEqual(hot['cv_vib_hartree_per_k'] / (Boltzmann / physical_constants['Hartree energy'][0]), 1, places=5)


@unittest.skipUnless(READY, 'Requires WSL PySCF and geomeTRIC')
class FrequencyPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.cfg = DFTConfig(atoms=parse_coordinates(WATER), xc='pbe', basis='sto-3g', engine='native',
                             grid_level=3, cube_points=12, threads=1)

    def tearDown(self):
        self.temp.cleanup()

    def optimized(self, cfg, name):
        from engines import calculate
        result = calculate(replace(cfg, calculation='optimize', opt_convergence='tight'), self.root / name,
                           lambda _: None)
        self.assertTrue(result['converged'])
        return DFTConfig.from_dict(result['final_geometry_config']).atoms

    def test_projection_and_units_match_pyscf(self):
        # Pure check of mass weighting, rigid-body removal and unit conversion on one Hessian.
        from pyscf import dft, gto
        from pyscf.hessian import thermo
        from engines.frequencies import harmonic_analysis, isotope_masses
        for atoms in (WATER, 'C 0 0 0; O 0 0 1.13', 'N 0 0 .1; H 0 .94 -.27; H .81 -.47 -.27; H -.81 -.47 -.27'):
            with self.subTest(atoms=atoms):
                mol = gto.M(atom=atoms, basis='sto-3g', verbose=0)
                mf = dft.RKS(mol)
                mf.xc = 'lda,vwn'
                mf.kernel()
                hessian = mf.Hessian().kernel()
                masses = isotope_masses([mol.atom_pure_symbol(i) for i in range(mol.natm)])
                reference = thermo.harmonic_analysis(mol, hessian, mass=masses)['freq_wavenumber']
                reference = np.where(np.iscomplex(reference), -np.abs(reference.imag), reference.real)
                flat = hessian.transpose(0, 2, 1, 3).reshape(3 * mol.natm, 3 * mol.natm)
                ours = harmonic_analysis(flat, mol.atom_coords(), masses)
                np.testing.assert_allclose(np.sort(ours['frequencies_cm1']), np.sort(reference), atol=1e-4)
                self.assertEqual(ours['linear'], mol.natm == 2)

    def test_water_matches_analytic_hessian(self):
        from engines import calculate
        from engines.frequencies import CM1_TO_HARTREE, isotope_masses
        from engines.pyscf_engine import configure_grids
        from pyscf import dft, gto
        from pyscf.hessian import thermo
        atoms = self.optimized(self.cfg, 'opt')
        result = calculate(replace(self.cfg, atoms=atoms, calculation='frequencies'), self.root / 'freq', lambda _: None)
        block = result['frequencies']
        self.assertEqual(result['status'], 'converged')
        self.assertEqual(block['vibrational_mode_count'], 3)
        self.assertEqual(block['imaginary_frequencies_cm1'], [])
        self.assertEqual(block['stationary_point_character'], 'minimum')
        self.assertEqual(block['gradient_calculations'], 19)
        self.assertEqual(len(list((self.root / 'freq' / 'displacements').iterdir())), 19)
        mol = gto.M(atom=[(a.element, a.position) for a in atoms], basis='sto-3g', verbose=0)
        mf = dft.RKS(mol)
        mf.xc, mf.small_rho_cutoff, mf.conv_tol = 'pbe', 0, 1e-11
        configure_grids(mf.grids, self.cfg)
        mf.kernel()
        reference = thermo.harmonic_analysis(mol, mf.Hessian().kernel(),
                                             mass=isotope_masses(['O', 'H', 'H']))['freq_wavenumber'].real
        np.testing.assert_allclose(block['frequencies_cm1'], np.sort(reference), atol=2.0)
        thermo_block = block['thermochemistry']
        self.assertAlmostEqual(thermo_block['zpe_hartree'], .5 * sum(block['frequencies_cm1']) * CM1_TO_HARTREE, places=12)
        modes = np.asarray(block['normal_modes_cartesian'])
        self.assertEqual(modes.shape, (3, 3, 3))
        np.testing.assert_allclose(np.linalg.norm(modes.reshape(3, -1), axis=1), 1, atol=1e-12)

    def test_h2_matches_energy_curvature(self):
        from engines import calculate
        from engines.frequencies import AU_TO_CM1
        cfg = replace(self.cfg, atoms=parse_coordinates('H 0 0 0; H 0 0 .75'))
        atoms = self.optimized(cfg, 'h2opt')
        result = calculate(replace(cfg, atoms=atoms, calculation='frequencies'), self.root / 'h2', lambda _: None)
        self.assertTrue(result['frequencies']['linear'])
        self.assertEqual(len(result['frequencies']['frequencies_cm1']), 1)
        r0 = np.linalg.norm(np.subtract(atoms[1].position, atoms[0].position))
        step, energies = 0.005, []
        for k in (-1, 0, 1):
            run = calculate(replace(cfg, atoms=parse_coordinates(f'H 0 0 0; H 0 0 {r0 + k * step * 0.529177210903}')),
                            self.root / f'e{k}', lambda _: None)
            energies.append(run['total_energy_hartree'])
        curvature = (energies[0] - 2 * energies[1] + energies[2]) / step**2
        expected = np.sqrt(curvature / (1.007825 / 2)) * AU_TO_CM1
        self.assertAlmostEqual(result['frequencies']['frequencies_cm1'][0], expected, delta=1.0)

    def test_ideal_gas_thermo_matches_pyscf(self):
        # Same Hessian, masses (PySCF thermo uses isotope-averaged masses) and frequencies:
        # entropies differ only by the CODATA constant versions (~4e-7 relative).
        from pyscf import dft, gto
        from pyscf.hessian import thermo
        from engines.frequencies import harmonic_analysis, ideal_gas_thermo, symmetry_number
        for atoms, spin, rotor in ((WATER, 0, 'nonlinear'), ('O 0 0 0; O 0 0 1.21', 2, 'linear')):
            with self.subTest(atoms=atoms):
                mol = gto.M(atom=atoms, basis='sto-3g', spin=spin, verbose=0)
                mf = dft.UKS(mol) if spin else dft.RKS(mol)
                mf.xc = 'lda,vwn'
                mf.kernel()
                hessian = mf.Hessian().kernel()
                masses = mol.atom_mass_list(isotope_avg=True)
                reference = thermo.thermo(mf, thermo.harmonic_analysis(mol, hessian, mass=masses)['freq_au'],
                                          298.15, 101325)
                flat = hessian.transpose(0, 2, 1, 3).reshape(3 * mol.natm, 3 * mol.natm)
                frequencies = harmonic_analysis(flat, mol.atom_coords(), masses)['frequencies_cm1']
                sigma = symmetry_number([mol.atom_pure_symbol(i) for i in range(mol.natm)], mol.atom_coords())
                ours = ideal_gas_thermo(mf.e_tot, frequencies, masses, mol.atom_coords(), spin + 1, sigma,
                                        298.15, 101325)
                self.assertEqual((ours['rotor'], sigma), (rotor, reference['sym_number'][0]))
                terms = ours['terms_hartree']
                for name, value in (('S_trans', terms['translational']['entropy']),
                                    ('S_rot', terms['rotational']['entropy']),
                                    ('S_elec', terms['electronic']['entropy'])):
                    self.assertAlmostEqual(value, reference[name][0], delta=1e-5 * max(abs(reference[name][0]), 1e-12))
                self.assertAlmostEqual(ours['enthalpy_hartree'], reference['H_tot'][0], delta=1e-7)
                self.assertAlmostEqual(ours['gibbs_free_energy_hartree'], reference['G_tot'][0], delta=1e-7)

    def test_gibbs_in_reaction_analysis_requires_same_conditions(self):
        from engines import calculate
        from energy_analysis import reaction_energy
        h2 = replace(self.cfg, atoms=parse_coordinates('H 0 0 0; H 0 0 .74'), calculation='frequencies', grid_level=2)
        h = replace(h2, atoms=parse_coordinates('H 0 0 0'), spin=1)
        mol_run = calculate(h2, self.root / 'h2', lambda _: None)
        atom_run = calculate(h, self.root / 'h', lambda _: None)
        for field in ('zpe_corrected_energy_hartree', 'enthalpy_hartree', 'gibbs_free_energy_hartree'):
            self.assertIn(field, mol_run)
        report = reaction_energy([(1, self.root / 'h2'), (-2, self.root / 'h')], 'gibbs_free_energy_hartree')
        expected = mol_run['gibbs_free_energy_hartree'] - 2 * atom_run['gibbs_free_energy_hartree']
        self.assertAlmostEqual(report['energy_hartree'], expected, places=12)
        self.assertEqual(report['conditions']['thermo_temperature_k'], 298.15)
        self.assertEqual(atom_run['frequencies']['ideal_gas_thermochemistry']['rotor'], 'atom')
        calculate(replace(h, thermo_temperature_k=500), self.root / 'h500', lambda _: None)
        with self.assertRaisesRegex(ValueError, 'thermo_temperature_k'):
            reaction_energy([(1, self.root / 'h2'), (-2, self.root / 'h500')], 'gibbs_free_energy_hartree')
        # The electronic energy of the same runs stays comparable regardless of temperature.
        reaction_energy([(1, self.root / 'h2'), (-2, self.root / 'h500')])

    def test_planar_ammonia_is_a_saddle_and_failures_stop(self):
        from engines import calculate
        planar = replace(self.cfg, atoms=parse_coordinates('N 0 0 0; H 0 1.01 0; H .8747 -.505 0; H -.8747 -.505 0'),
                         calculation='frequencies', grid_level=2)
        result = calculate(planar, self.root / 'nh3', lambda _: None)
        block = result['frequencies']
        self.assertGreaterEqual(len([f for f in block['imaginary_frequencies_cm1'] if f < -20]), 1)
        self.assertNotEqual(block['stationary_point_character'], 'minimum')
        umbrella = np.asarray(block['normal_modes_cartesian'][0])
        self.assertGreater(np.linalg.norm(umbrella[:, 2]), .99)  # Out of the molecular plane.
        self.assertTrue(any('허수 진동' in w for w in result['warnings']))
        with self.assertRaisesRegex(RuntimeError, '변위 계산의 SCF'):
            calculate(replace(planar, max_cycle=1), self.root / 'failed', lambda _: None)


if __name__ == '__main__':
    unittest.main()
