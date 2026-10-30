"""Analytic response vs Fock/gradient differences and independent PySCF response."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from molecular import DFTConfig, parse_coordinates


@unittest.skipUnless(importlib.util.find_spec('pyscf'), 'Requires WSL PySCF')
class NativeResponseTests(unittest.TestCase):
    def test_rks_uks_lda_gga_hybrid_response(self):
        from pyscf import gto, dft, lib
        from engines.native_scf import NativeKS
        from engines.native_response import FockResponse
        from engines.native_stability import OrbitalRotations
        from engines.pyscf_engine import configure_grids
        old_threads = lib.num_threads()
        lib.num_threads(1)
        self.addCleanup(lib.num_threads, old_threads)
        rng = np.random.default_rng(23)
        with tempfile.TemporaryDirectory() as temp:
            for atoms, spin in [('O 0 0 0; H 0 -.757 .587; H 0 .757 .587', 0),
                                ('N 0 0 0; H 0 -.8 .6; H 0 .8 .6', 1)]:
                for xc in ('lda,vwn', 'pbe', 'pbe0', 'b3lyp'):
                    with self.subTest(spin=spin, xc=xc):
                        cfg = DFTConfig(atoms=parse_coordinates(atoms), spin=spin, basis='sto-3g',
                                        xc=xc, engine='native', grid_level=1, threads=1)
                        mol = gto.M(atom=atoms, spin=spin, basis=cfg.basis, verbose=0)
                        mf = NativeKS(mol, cfg, Path(temp))
                        configure_grids(mf.grids, cfg)
                        mf.grids.build()
                        mf.mo_energy, mf.mo_coeff = mf.diagonalize(
                            np.stack((mf.hcore, mf.hcore)) if spin else mf.hcore)
                        initial = OrbitalRotations(mf)
                        rotated = initial.rotated(.1 * rng.normal(size=initial.dimension))
                        mf.mo_coeff = np.asarray(rotated) if spin else rotated[0]
                        dm = mf.density(mf.mo_coeff)
                        rotations = OrbitalRotations(mf)
                        direction = rng.normal(size=rotations.dimension)
                        direction /= np.linalg.norm(direction)
                        h = 1e-5
                        plus = rotations.density(rotations.rotated(h * direction))
                        minus = rotations.density(rotations.rotated(-h * direction))
                        delta = (plus - minus) / (2 * h)
                        response = FockResponse(mf, dm)
                        actual = response(delta)
                        # Stay on physical, positive densities; arbitrary additive D
                        # perturbations can cross Libxc's low-density clipping threshold.
                        finite = (mf.evaluate(plus)[1] - mf.evaluate(minus)[1]) / (2 * h)
                        np.testing.assert_allclose(actual, finite, atol=2e-7, rtol=2e-5)
                        ref = dft.UKS(mol) if spin else dft.RKS(mol)
                        ref.xc, ref.grids = xc, mf.grids
                        # PySCF builds the entire response only in this independent test.
                        reference = ref.gen_response(mf.mo_coeff, mf.mo_occ, hermi=1)(delta)
                        np.testing.assert_allclose(actual, reference, atol=2e-11, rtol=2e-10)
                        response.cache = None  # Streaming fallback must use the same quadrature.
                        np.testing.assert_allclose(response(delta), actual, atol=1e-13)
                        x, y = rng.normal(size=(2, rotations.dimension))
                        x, y = x / np.linalg.norm(x), y / np.linalg.norm(y)
                        hx, hy = rotations.hessian_vector(x), rotations.hessian_vector(y)
                        np.testing.assert_allclose(hx, rotations.hessian_vector_fd(x, step=1e-5),
                                                   atol=2e-7, rtol=2e-6)
                        self.assertAlmostEqual(x @ hy, y @ hx, delta=2e-10)
                        # Cached Hessian applications must not create new Fock matrices.
                        with patch.object(mf, 'evaluate', side_effect=AssertionError('Fock rebuild')):
                            np.testing.assert_allclose(rotations.hessian_vector(2 * x), 2 * hx, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
