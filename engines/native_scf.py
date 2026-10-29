"""Inspectable molecular KS solver; no PySCF SCF/DFT driver is used.

External primitives: PySCF/libcint integrals and AO values, atom-centred
quadrature, Libxc local functional derivatives, NumPy/SciPy linear algebra.
Density, gradient, XC matrix assembly, J/K contractions, occupations, Pulay
DIIS, SCF stopping and energy accounting are implemented in this module.
"""
import numpy as np
from scipy.linalg import eigh, expm
from scf_control import check_cancelled


class PulayDIIS:
    def __init__(self, size):
        self.size, self.focks, self.errors = size, [], []

    def update(self, fock, error):
        self.focks.append(fock.copy())
        self.errors.append(error.ravel().copy())
        self.focks, self.errors = self.focks[-self.size:], self.errors[-self.size:]
        n = len(self.errors)
        if n < 2:
            return fock
        gram = np.array([[np.dot(a, b) for b in self.errors] for a in self.errors])
        scale = max(float(np.max(np.abs(gram))), 1e-30)
        matrix = np.empty((n + 1, n + 1))
        matrix[:n, :n] = gram / scale
        matrix[:n, n] = matrix[n, :n] = -1
        matrix[n, n] = 0
        rhs = np.zeros(n + 1)
        rhs[-1] = -1
        try:
            weights = np.linalg.solve(matrix, rhs)[:n]
            if not np.isfinite(weights).all() or np.max(np.abs(weights)) > 30:
                raise np.linalg.LinAlgError('Ill-conditioned Pulay history')
            return np.einsum('i,i...->...', weights, np.asarray(self.focks))
        except np.linalg.LinAlgError:
            self.focks, self.errors = self.focks[-1:], self.errors[-1:]
            return fock


class NativeKS:
    def __init__(self, mol, config, directory):
        from pyscf.dft import gen_grid, numint, libxc
        self.mol, self.config, self.directory = mol, config, directory
        self.xc = config.xc
        self.grids = gen_grid.Grids(mol)
        # Used only by the common exporter for an independent grid-density check.
        self._numint = numint.NumInt()
        self.xctype = libxc.xc_type(self.xc)
        if self.xctype not in ('LDA', 'GGA') or libxc.rsh_coeff(self.xc)[0] != 0:
            raise ValueError('자체 엔진은 LDA/GGA와 지정된 global hybrid만 지원합니다.')
        self.hybrid = float(libxc.hybrid_coeff(self.xc))
        self.overlap = mol.intor_symmetric('int1e_ovlp')
        eigenvalues, vectors = eigh(self.overlap)
        if eigenvalues.min() < 1e-9:
            raise ValueError('겹침행렬이 거의 선형종속입니다. 더 안정적인 기저를 선택하세요.')
        self.orthogonalizer = vectors / np.sqrt(eigenvalues)
        self.kinetic = mol.intor_symmetric('int1e_kin')
        self.attraction = mol.intor_symmetric('int1e_nuc')
        self.hcore = self.kinetic + self.attraction
        n = mol.nao_nr()
        # Reserve memory for contractions, AO chunks, and eigensolvers too.
        if n > 160 or 3 * n**4 * 8 + 64 * 1024**2 > config.memory_mb * 1024**2:
            raise ValueError('자체 엔진의 dense ERI는 최대160 AO이며 메모리 여유가 필요합니다. '
                             '작은 기저/더 큰 메모리 목표를 사용하거나 PySCF 엔진을 선택하세요.')
        check_cancelled(directory)
        self.eri = mol.intor('int2e', aosym='s1')
        self.unrestricted = config.resolved_method == 'uks'
        self.mo_occ = np.zeros((2, n) if self.unrestricted else n)
        if self.unrestricted:
            for spin, electrons in enumerate(mol.nelec):
                self.mo_occ[spin, :electrons] = 1
        else:
            self.mo_occ[:mol.nelectron // 2] = 2
        self.converged = False
        self.callback = lambda env: None
        self.max_cycle = config.max_cycle
        self.conv_tol = config.conv_tol
        self.scf_summary = {}

    def diagonalize(self, fock):
        matrices = fock if fock.ndim == 3 else [fock]
        solutions = [eigh(self.orthogonalizer.T @ f @ self.orthogonalizer) for f in matrices]
        energies = np.asarray([e for e, _ in solutions])
        coeff = np.asarray([self.orthogonalizer @ c for _, c in solutions])
        return (energies, coeff) if fock.ndim == 3 else (energies[0], coeff[0])

    def density(self, coeff):
        return (coeff * self.mo_occ[..., None, :]) @ coeff.swapaxes(-1, -2)

    def xc_matrix(self, dm):
        from pyscf.dft import numint, libxc
        matrices = dm if self.unrestricted else dm[None]
        potentials = np.zeros_like(matrices)
        exc = 0.0
        gga = self.xctype == 'GGA'
        for start in range(0, len(self.grids.weights), 2048):
            check_cancelled(self.directory)
            coords = self.grids.coords[start:start + 2048]
            weights = self.grids.weights[start:start + 2048]
            ao = numint.eval_ao(self.mol, coords, deriv=int(gga))
            values = ao[0] if gga else ao
            rho = []
            for matrix in matrices:
                product = values @ matrix
                density = np.einsum('gi,gi->g', product, values)
                if gga:
                    gradient = 2 * np.einsum('gi,xgi->xg', product, ao[1:4])
                    rho.append(np.vstack((density, gradient)))
                else:
                    rho.append(density)
            rho = np.asarray(rho)
            eps, vxc, _, _ = libxc.eval_xc(self.xc, rho if self.unrestricted else rho[0],
                                          spin=int(self.unrestricted), deriv=1)
            scalar_density = rho[:, 0] if gga else rho
            exc += float(np.dot(weights * scalar_density.sum(axis=0), eps))
            vrho = vxc[0].T if self.unrestricted else vxc[0][None]
            for spin in range(len(matrices)):
                potentials[spin] += values.T @ ((weights * vrho[spin])[:, None] * values)
                if gga:
                    if self.unrestricted:
                        vsigma = vxc[1].T
                        derivative = 2 * vsigma[0 if spin == 0 else 2] * rho[spin, 1:4]
                        derivative += vsigma[1] * rho[1 - spin, 1:4]
                    else:
                        derivative = 2 * vxc[1] * rho[0, 1:4]
                    for axis in range(3):
                        term = ao[axis + 1].T @ ((weights * derivative[axis])[:, None] * values)
                        potentials[spin] += term + term.T
        return exc, potentials if self.unrestricted else potentials[0]

    def evaluate(self, dm):
        total_dm = dm.sum(axis=0) if self.unrestricted else dm
        coulomb = np.einsum('pqrs,rs->pq', self.eri, total_dm, optimize=True)
        exc, vxc = self.xc_matrix(dm)
        exchange = np.einsum('prqs,...rs->...pq', self.eri, dm, optimize=True) if self.hybrid else np.zeros_like(dm)
        factor = 1 if self.unrestricted else .5
        fock = self.hcore + coulomb + vxc - factor * self.hybrid * exchange
        exx = -factor * .5 * self.hybrid * float(np.sum(dm * exchange.swapaxes(-1, -2)))
        components = {'kinetic': float(np.einsum('ij,ji->', total_dm, self.kinetic)),
                      'nuclear_attraction': float(np.einsum('ij,ji->', total_dm, self.attraction)),
                      'coulomb': .5 * float(np.einsum('ij,ji->', total_dm, coulomb)),
                      'xc_semilocal': exc, 'exact_exchange': exx,
                      'nuclear_repulsion': float(self.mol.energy_nuc())}
        return sum(components.values()), fock, components

    def residual(self, fock, dm):
        error = fock @ dm @ self.overlap - self.overlap @ dm @ fock
        return self.orthogonalizer.T @ error @ self.orthogonalizer

    def kernel(self, dm0=None):
        self.grids.build()
        h = np.stack((self.hcore, self.hcore)) if self.unrestricted else self.hcore
        self.mo_energy, self.mo_coeff = self.diagonalize(h)
        if self.config.initial_guess == 'perturbed' and dm0 is None:
            rng = np.random.default_rng(1729)
            random = rng.normal(size=self.hcore.shape)
            rotation = expm(.08 * (random - random.T))
            self.mo_coeff = self.mo_coeff @ rotation
        dm = self.density(self.mo_coeff) if dm0 is None else dm0.copy()
        energy, fock, components = self.evaluate(dm)
        diis = PulayDIIS(self.config.native_diis_space)
        for cycle in range(self.max_cycle):
            check_cancelled(self.directory)
            extrapolated = diis.update(fock, self.residual(fock, dm))
            energies, coeff = self.diagonalize(extrapolated)
            pure = self.density(coeff)
            trial = (1 - self.config.native_damping) * pure + self.config.native_damping * dm
            next_energy, next_fock, next_components = self.evaluate(trial)
            de = next_energy - energy
            ddm = float(np.linalg.norm(trial - dm))
            residual = float(np.linalg.norm(self.residual(next_fock, trial)))
            ok = abs(de) < self.conv_tol and ddm < self.config.native_density_tol and residual < self.config.native_residual_tol
            if ok:
                # Check canonical, integer-occupation output against the same
                # criteria. No relaxed "extra cycle" success thresholds.
                energies, coeff = self.diagonalize(next_fock)
                pure = self.density(coeff)
                canonical_energy, canonical_fock, canonical_parts = self.evaluate(pure)
                ok = (abs(canonical_energy - next_energy) < self.conv_tol
                      and abs(canonical_energy - energy) < self.conv_tol
                      and np.linalg.norm(pure - trial) < self.config.native_density_tol
                      and np.linalg.norm(pure - dm) < self.config.native_density_tol
                      and np.linalg.norm(self.residual(canonical_fock, pure)) < self.config.native_residual_tol)
                if ok:
                    trial, next_energy, next_fock, next_components = pure, canonical_energy, canonical_fock, canonical_parts
                    de = next_energy - energy
                    ddm = float(np.linalg.norm(trial - dm))
                    residual = float(np.linalg.norm(self.residual(next_fock, trial)))
            self.mo_energy, self.mo_coeff = energies, coeff
            self._dm = trial
            self.callback({'cycle': cycle, 'e_tot': next_energy, 'norm_gorb': residual,
                           'norm_ddm': ddm, 'delta_energy': de, 'dm': trial,
                           'components': next_components})
            dm, energy, fock, components = trial, next_energy, next_fock, next_components
            if ok:
                self.converged = True
                break
        if not self.converged:
            # Save a self-consistent representation of the last orbital trial,
            # but never promote it to convergence or claim another iteration.
            self._dm = self.density(self.mo_coeff)
            energy, fock, components = self.evaluate(self._dm)
        self.components = components
        self.final_residual = float(np.linalg.norm(self.residual(fock, self._dm)))
        self.scf_summary = {'e1': components['kinetic'] + components['nuclear_attraction'],
                            'coul': components['coulomb'],
                            'exc': components['xc_semilocal'] + components['exact_exchange'],
                            'nuc': components['nuclear_repulsion']}
        return energy

    def get_hcore(self):
        return self.hcore

    def make_rdm1(self):
        return self._dm.copy()

    def spin_square(self):
        if not self.unrestricted:
            return 0., 1.
        alpha, beta = self.mol.nelec
        cross = self.mo_coeff[0, :, :alpha].T @ self.overlap @ self.mo_coeff[1, :, :beta]
        sz = .5 * (alpha - beta)
        s2 = float(sz * (sz + 1) + beta - np.sum(cross**2))
        return s2, float(np.sqrt(max(0., 1 + 4 * s2)))
