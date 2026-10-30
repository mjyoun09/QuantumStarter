"""Real internal SCF stability of native RKS/UKS solutions, and following instabilities.

Orbital rotations C(k) = C exp(K), K[a,i] = k_ai = -K[i,a] (occupied i, virtual a,
per spin). The energy gradient with respect to k at the current orbitals is
g_ai = 2 n F_ai with n = 2 (RKS) or 1 (UKS) and F in the MO basis. Hessian-vector
products use a native analytic Fock response with the Libxc XC kernel. Central
differences are retained as an independent validation path.
The lowest Hessian eigenvalue comes from Lanczos (scipy eigsh). Only real
internal rotations are tested: RKS -> UKS (triplet) instabilities are not.
"""
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import LinearOperator, eigsh

from scf_control import check_cancelled
from .native_response import FockResponse


class OrbitalRotations:
    def __init__(self, mf):
        self.mf = mf
        coeff = mf.mo_coeff if mf.unrestricted else mf.mo_coeff[None]
        occ = mf.mo_occ if mf.unrestricted else mf.mo_occ[None]
        self.coeff = [np.asarray(c) for c in coeff]
        self.occupied = [np.flatnonzero(o > 1e-8) for o in occ]
        self.virtual = [np.flatnonzero(o <= 1e-8) for o in occ]
        self.occupation = 1.0 if mf.unrestricted else 2.0
        self.sizes = [len(o) * len(v) for o, v in zip(self.occupied, self.virtual)]
        self._response = None
        self._focks = None

    @property
    def dimension(self):
        return sum(self.sizes)

    def rotated(self, kappa):
        coefficients, start = [], 0
        for c, occ, vir, size in zip(self.coeff, self.occupied, self.virtual, self.sizes):
            block = np.zeros((c.shape[1], c.shape[1]))
            k = kappa[start:start + size].reshape(len(vir), len(occ))
            block[np.ix_(vir, occ)] = k
            block[np.ix_(occ, vir)] = -k.T
            coefficients.append(c @ expm(block))
            start += size
        return coefficients

    def density(self, coefficients):
        densities = [self.occupation * c[:, occ] @ c[:, occ].T for c, occ in zip(coefficients, self.occupied)]
        return np.asarray(densities) if self.mf.unrestricted else densities[0]

    def energy_and_gradient(self, kappa):
        check_cancelled(self.mf.directory)
        coefficients = self.rotated(kappa)
        energy, fock, _ = self.mf.evaluate(self.density(coefficients))
        focks = fock if self.mf.unrestricted else fock[None]
        parts = [2 * self.occupation * (c.T @ f @ c)[np.ix_(vir, occ)].ravel()
                 for c, f, occ, vir in zip(coefficients, focks, self.occupied, self.virtual)]
        return energy, np.concatenate(parts)

    def _prepare_response(self):
        if self._response is None:
            density = self.density(self.coeff)
            _, fock, _ = self.mf.evaluate(density)
            self._focks = fock if self.mf.unrestricted else fock[None]
            self._response = FockResponse(self.mf, density)

    def diagonal(self):
        """Positive orbital-gap preconditioner, including spin occupation factors."""
        self._prepare_response()
        parts = []
        for c, f, occ, vir in zip(self.coeff, self._focks, self.occupied, self.virtual):
            energies = np.diag(c.T @ f @ c)
            parts.append((2 * self.occupation * (energies[vir, None] - energies[None, occ])).ravel())
        return np.maximum(np.abs(np.concatenate(parts)), 1e-3)

    def hessian_vector(self, vector):
        """d(2 n F_ai)/dk, with dF = J[dD] + f_xc[dD] - hybrid*K[dD]."""
        if not np.any(vector):
            return np.zeros_like(vector)
        self._prepare_response()
        delta_coeff, delta_dm, start = [], [], 0
        for c, occ, vir, size in zip(self.coeff, self.occupied, self.virtual, self.sizes):
            k = vector[start:start + size].reshape(len(vir), len(occ))
            dc = np.zeros_like(c)
            dc[:, occ] = c[:, vir] @ k
            dc[:, vir] = -c[:, occ] @ k.T
            term = self.occupation * dc[:, occ] @ c[:, occ].T
            delta_dm.append(term + term.T)
            delta_coeff.append(dc)
            start += size
        response = self._response(np.asarray(delta_dm) if self.mf.unrestricted else delta_dm[0])
        responses = response if self.mf.unrestricted else response[None]
        parts = []
        for c, dc, f, df, occ, vir in zip(self.coeff, delta_coeff, self._focks, responses,
                                         self.occupied, self.virtual):
            derivative = dc.T @ f @ c + c.T @ f @ dc + c.T @ df @ c
            parts.append((2 * self.occupation * derivative[np.ix_(vir, occ)]).ravel())
        return np.concatenate(parts)

    def hessian_vector_fd(self, vector, step=1e-4):
        norm = np.linalg.norm(vector)
        if norm == 0:
            return np.zeros_like(vector)
        direction = vector / norm
        plus = self.energy_and_gradient(step * direction)[1]
        minus = self.energy_and_gradient(-step * direction)[1]
        return norm * (plus - minus) / (2 * step)


def lowest_mode(rotations, tolerance=1e-6):
    if rotations.dimension == 0:
        return 0.0, np.zeros(0)
    if rotations.dimension == 1:
        return float(rotations.hessian_vector(np.ones(1))[0]), np.ones(1)
    operator = LinearOperator((rotations.dimension,) * 2, matvec=rotations.hessian_vector, dtype=float)
    start = np.random.default_rng(7).normal(size=rotations.dimension)
    values, vectors = eigsh(operator, k=1, which='SA', tol=tolerance, v0=start,
                            ncv=min(rotations.dimension, 20), maxiter=200)
    return float(values[0]), vectors[:, 0]


def stability_loop(mf, config, progress, max_rounds=3, threshold=-1e-4):
    """Check the converged solution; with 'follow', rotate downhill and re-run the SCF."""
    rounds, energy = [], float(mf.e_tot)
    for round_index in range(max_rounds):
        rotations = OrbitalRotations(mf)
        if rotations.dimension == 0:
            rounds.append({'round': round_index + 1, 'lowest_eigenvalue_hartree': None, 'energy_hartree': energy})
            break
        value, vector = lowest_mode(rotations)
        row = {'round': round_index + 1, 'lowest_eigenvalue_hartree': value, 'energy_hartree': energy,
               'dimension': rotations.dimension, 'stable': value > threshold}
        rounds.append(row)
        progress(f'SCF 안정성 {round_index + 1}회 · 최저 궤도 Hessian 고유값 {value:+.3e} Ha · '
                 f'{"안정" if row["stable"] else "불안정(더 낮은 해 존재)"}')
        if row['stable'] or config.native_stability != 'follow':
            break
        # Line search along the unstable mode, then a fresh SCF from the rotated orbitals.
        trials = [(rotations.energy_and_gradient(t * vector)[0], t) for t in (.05, .1, .2, .4, .8, -.1, -.4)]
        best, step = min(trials)
        row['line_search'] = {'step_radian': step, 'energy_hartree': best}
        if best >= energy:
            row['note'] = 'line search found no lower energy along the unstable mode'
            break
        guess = rotations.density(rotations.rotated(step * vector))
        energy = float(mf.kernel(dm0=guess))
        if not mf.converged:
            row['note'] = 'SCF from the rotated orbitals did not converge'
            break
    final = rounds[-1]
    return energy, {'mode': config.native_stability, 'kind': 'real internal (RKS->RKS / UKS->UKS)',
                    'threshold_hartree': threshold, 'rounds': rounds,
                    'stable': bool(final.get('stable', True)) and final.get('note') is None
                    and final['energy_hartree'] == energy,
                    'not_tested': 'RKS->UKS (triplet) and complex instabilities'}
