"""Native trust-region Newton SCF: same minimum as DIIS, escape from an unstable (saddle) solution."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from molecular import DFTConfig, parse_coordinates


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class NativeNewtonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_water_matches_diis(self):
        from engines import calculate
        cfg = DFTConfig(atoms=parse_coordinates('O 0 0 0; H 0 -.757 .587; H 0 .757 .587'), xc='pbe', basis='sto-3g',
                        engine='native', scf_solver='newton', grid_level=2, cube_points=12, threads=1)
        newton = calculate(cfg, self.root / 'newton', lambda _: None)
        diis = calculate(replace(cfg, scf_solver='diis'), self.root / 'diis', lambda _: None)
        self.assertTrue(newton['converged'])
        self.assertAlmostEqual(newton['total_energy_hartree'], diis['total_energy_hartree'], delta=1e-8)
        self.assertLess(newton['scf_criteria']['final_residual'], cfg.native_residual_tol)
        self.assertLessEqual(newton['native_newton']['aufbau_violation_hartree'], 0)
        self.assertAlmostEqual(newton['electron_count_from_density'], 10, places=8)

    def test_escapes_unstable_solution(self):
        # DIIS converges OH (coarse grid) to a state that the stability analysis shows to be a
        # saddle point. Started there, Newton must leave it and stop only at a local minimum.
        from pyscf import gto
        from engines.native_scf import NativeKS
        from engines.native_stability import OrbitalRotations, lowest_mode
        from engines.pyscf_engine import configure_grids
        cfg = DFTConfig(atoms=parse_coordinates('O 0 0 0; H 0.1 0 .97'), spin=1, xc='pbe', basis='sto-3g',
                        engine='native', grid_level=1, conv_tol=1e-12, native_density_tol=1e-10,
                        native_residual_tol=1e-9, max_cycle=300)
        mol = gto.M(atom='O 0 0 0; H 0.1 0 .97', basis='sto-3g', spin=1, verbose=0)
        diis = NativeKS(mol, cfg, self.root)
        configure_grids(diis.grids, cfg)
        saddle = diis.kernel()
        self.assertLess(lowest_mode(OrbitalRotations(diis))[0], -1e-4)
        newton = NativeKS(mol, replace(cfg, scf_solver='newton', max_cycle=100), self.root)
        configure_grids(newton.grids, cfg)
        energy = newton.kernel(dm0=diis.make_rdm1())
        self.assertTrue(newton.converged)
        self.assertLess(energy, saddle - 1e-7)
        self.assertTrue(any(row['lowest_eigenvalue_hartree'] < -1e-4 for row in newton.newton_escapes))
        self.assertGreater(lowest_mode(OrbitalRotations(newton))[0], -1e-4)
        self.assertAlmostEqual(np.trace(newton.make_rdm1().sum(axis=0) @ newton.overlap), 9, places=8)


if __name__ == '__main__':
    unittest.main()
