"""Native fixed-grid Fock response dF/dD for LDA/GGA/global hybrids.

Only AO values and Libxc local derivatives (in Cartesian density/gradient
variables) come from PySCF. Density response, quadrature matrices and J/K are
assembled here; no PySCF response/SCF driver is used. Cache within a memory
budget, otherwise stream chunks on each application.
"""
import numpy as np

from scf_control import check_cancelled


def grid_density(ao, matrices):
    """Spin density and its Cartesian derivatives for symmetric real D."""
    values = ao[0]
    rho = []
    for matrix in matrices:
        product = values @ matrix
        density = np.einsum('gi,gi->g', product, values)
        gradient = 2 * np.einsum('gi,xgi->xg', product, ao[1:])
        rho.append(np.vstack((density, gradient)))
    return np.asarray(rho)


class FockResponse:
    def __init__(self, mf, density):
        self.mf = mf
        self.density = density.copy()
        nspin, nvar = (2 if mf.unrestricted else 1), (4 if mf.xctype == 'GGA' else 1)
        estimate = len(mf.grids.weights) * (nvar * mf.mol.nao_nr() + (nspin * nvar)**2) * 8
        available = max(0, mf.config.memory_mb * 1024**2 - 3 * mf.eri.nbytes - 64 * 1024**2)
        self.cache = list(self._chunks()) if estimate < available / 2 else None

    def _chunks(self):
        from pyscf.dft import numint, libxc
        mf = self.mf
        matrices = self.density if mf.unrestricted else self.density[None]
        gga = mf.xctype == 'GGA'
        for start in range(0, len(mf.grids.weights), 2048):
            check_cancelled(mf.directory)
            coords = mf.grids.coords[start:start + 2048]
            ao = numint.eval_ao(mf.mol, coords, deriv=int(gga))
            ao = ao if gga else ao[None]
            rho = grid_density(ao, matrices)
            # libxc.eval_xc_eff returns the requested derivative tensor itself,
            # unlike NumInt.eval_xc_eff, which returns a tuple of derivative orders.
            fxc = libxc.eval_xc_eff(mf.xc, rho if mf.unrestricted else rho[0],
                                    deriv=2, spin=int(mf.unrestricted))
            if not mf.unrestricted:
                fxc = fxc[None, :, None, :, :]
            yield ao, mf.grids.weights[start:start + 2048], fxc

    def __call__(self, delta_density, *, transition=False):
        mf = self.mf
        check_cancelled(mf.directory)
        matrices = delta_density if mf.unrestricted else delta_density[None]
        # A TDA transition density is oriented. Local density and gradient see its
        # symmetric part; exact exchange must contract the original matrix.
        local = (matrices + matrices.swapaxes(-1, -2)) / 2 if transition else matrices
        potentials = np.zeros_like(matrices)
        for ao, weights, fxc in self.cache if self.cache is not None else self._chunks():
            check_cancelled(mf.directory)
            delta_rho = grid_density(ao, local)
            delta_v = np.einsum('sxtyg,tyg->sxg', fxc, delta_rho, optimize=True)
            values = ao[0]
            for spin in range(len(matrices)):
                potentials[spin] += values.T @ ((weights * delta_v[spin, 0])[:, None] * values)
                for axis in range(1, len(ao)):
                    term = ao[axis].T @ ((weights * delta_v[spin, axis])[:, None] * values)
                    potentials[spin] += term + term.T
        total = matrices.sum(axis=0)
        coulomb = np.einsum('pqrs,rs->pq', mf.eri, total, optimize=True)
        potentials += coulomb
        if mf.hybrid:
            exchange = np.einsum('prqs,...rs->...pq', mf.eri, matrices, optimize=True)
            potentials -= (1 if mf.unrestricted else .5) * mf.hybrid * exchange
        return potentials if mf.unrestricted else potentials[0]
