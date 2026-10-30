"""Native SCF stability: orbital gradient/Hessian consistency, verdicts vs PySCF, following."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from molecular import DFTConfig, parse_coordinates

OH = 'O 0 0 0; H 0.1 0 .97'
WATER = 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587'


class StabilityInputTests(unittest.TestCase):
    def test_native_only(self):
        self.assertEqual(DFTConfig(engine='native', native_stability='follow').native_stability, 'follow')
        for kwargs in ({'native_stability': 'follow'}, {'engine': 'native', 'native_stability': 'yes'}):
            with self.subTest(**kwargs), self.assertRaises(ValueError):
                DFTConfig(**kwargs)


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class StabilityPhysicsTests(unittest.TestCase):
    def setUp(self):
        from pyscf import lib
        old_threads = lib.num_threads()
        lib.num_threads(1)
        self.addCleanup(lib.num_threads, old_threads)
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def config(self, atoms, spin, xc):
        return DFTConfig(atoms=parse_coordinates(atoms), spin=spin, xc=xc, basis='sto-3g', engine='native',
                         grid_level=1, conv_tol=1e-12, native_density_tol=1e-10, native_residual_tol=1e-9,
                         max_cycle=300, cube_points=12, threads=1)

    def solver(self, cfg, converge=True):
        from pyscf import gto
        from engines.native_scf import NativeKS
        from engines.pyscf_engine import configure_grids
        mol = gto.M(atom=[(a.element, a.position) for a in cfg.atoms], basis=cfg.basis, spin=cfg.spin, verbose=0)
        mf = NativeKS(mol, cfg, self.root)
        configure_grids(mf.grids, cfg)
        if converge:
            mf.kernel()
        else:  # Core-Hamiltonian orbitals: a deliberately non-stationary point.
            mf.grids.build()
            mf.mo_energy, mf.mo_coeff = mf.diagonalize(
                np.stack((mf.hcore, mf.hcore)) if mf.unrestricted else mf.hcore)
        return mf

    def test_gradient_and_hessian_match_energy_derivatives(self):
        from engines.native_stability import OrbitalRotations
        rng = np.random.default_rng(11)
        for atoms, spin, xc in ((WATER, 0, 'pbe0'), (OH, 1, 'pbe')):
            with self.subTest(atoms=atoms):
                cfg = self.config(atoms, spin, xc)
                rotations = OrbitalRotations(self.solver(cfg, converge=False))
                x = rng.normal(size=rotations.dimension)
                x /= np.linalg.norm(x)
                e0, g0 = rotations.energy_and_gradient(np.zeros(rotations.dimension))
                step = 1e-4
                plus = rotations.energy_and_gradient(step * x)[0]
                minus = rotations.energy_and_gradient(-step * x)[0]
                self.assertGreater(np.linalg.norm(g0), 1e-2)
                self.assertAlmostEqual((plus - minus) / (2 * step), g0 @ x, delta=1e-6 * max(1, abs(g0 @ x)))
                converged = OrbitalRotations(self.solver(cfg))
                e0 = converged.energy_and_gradient(np.zeros(converged.dimension))[0]
                plus = converged.energy_and_gradient(step * x)[0]
                minus = converged.energy_and_gradient(-step * x)[0]
                curvature = (plus - 2 * e0 + minus) / step**2
                self.assertAlmostEqual(x @ converged.hessian_vector(x), curvature, delta=1e-4 * abs(curvature))

    def test_verdicts_agree_with_pyscf_and_follow_lowers_energy(self):
        from engines import calculate
        from engines.pyscf_engine import configure_grids
        from pyscf import dft, gto
        for atoms, spin, xc in ((WATER, 0, 'pbe0'), (OH, 1, 'pbe')):
            with self.subTest(atoms=atoms):
                cfg = self.config(atoms, spin, xc)
                mf = self.solver(cfg)
                mol = mf.mol
                reference = dft.UKS(mol) if spin else dft.RKS(mol)
                reference.xc, reference.small_rho_cutoff, reference.max_cycle = xc, 0, 0
                configure_grids(reference.grids, cfg)
                reference.kernel(dm0=mf.make_rdm1())
                reference.mo_coeff, reference.mo_energy, reference.mo_occ = mf.mo_coeff, mf.mo_energy, mf.mo_occ
                status = (reference.stability(return_status=True) if spin else
                          reference.stability(internal=True, external=False, return_status=True))[2]
                checked = calculate(replace(cfg, native_stability='check'), self.root / f'check{spin}', lambda _: None)
                self.assertEqual(checked['scf_stability']['stable'], bool(status))
                if status:
                    self.assertGreater(checked['scf_stability']['rounds'][0]['lowest_eigenvalue_hartree'], 0)
                    continue
                self.assertTrue(any('SCF 안정성' in w for w in checked['warnings']))
                followed = calculate(replace(cfg, native_stability='follow'), self.root / f'follow{spin}', lambda _: None)
                block = followed['scf_stability']
                self.assertTrue(followed['converged'])
                self.assertTrue(block['stable'])
                self.assertGreater(len(block['rounds']), 1)
                self.assertLess(followed['total_energy_hartree'], checked['total_energy_hartree'] - 1e-8)
                self.assertGreater(block['rounds'][-1]['lowest_eigenvalue_hartree'], 0)


if __name__ == '__main__':
    unittest.main()
