"""Second-order (trust-region Newton) SCF for the native engine.

Minimizes E over real orbital rotations C exp(K) with the occupation pattern held
by orbital index, so the solution need not obey the aufbau rule (needed for open
d shells such as the Cr and Co atoms, where DIIS with aufbau occupations
oscillates). Steps solve (H) x = -g with Steihaug truncated conjugate gradients
inside a trust radius; g and Hessian-vector products come from
native_stability.OrbitalRotations (analytic gradient and H.v). An orbital-gap
preconditioner resolves the stiff core and soft valence directions together.
Convergence also requires the lowest real internal curvature >= -1e-4 Ha;
this numerical check does not establish a global ground state.
"""
import numpy as np
from scipy.linalg import eigh

from scf_control import check_cancelled
from .native_stability import OrbitalRotations, lowest_mode


def steihaug(hessian_vector, gradient, radius, tolerance, max_iterations=100, diagonal=None):
    """Approximate argmin g.x + x.Hx/2 subject to |x| <= radius."""
    x = np.zeros_like(gradient)
    residual = -gradient
    diagonal = np.ones_like(gradient) if diagonal is None else diagonal
    z = residual / diagonal
    direction = z.copy()
    rz = residual @ z
    if np.linalg.norm(residual) <= tolerance:
        return x

    def to_boundary(x, d):
        a, b, c = d @ d, 2 * x @ d, x @ x - radius**2
        return x + (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a) * d

    for _ in range(max_iterations):
        hd = hessian_vector(direction)
        curvature = direction @ hd
        if curvature <= 0:  # Negative curvature: go to the trust boundary along it.
            return to_boundary(x, direction)
        alpha = rz / curvature
        trial = x + alpha * direction
        if np.linalg.norm(trial) >= radius:
            return to_boundary(x, direction)
        new_residual = residual - alpha * hd
        if np.linalg.norm(new_residual) < tolerance:
            return trial
        z = new_residual / diagonal
        next_rz = new_residual @ z
        direction = z + next_rz / rz * direction
        rz = next_rz
        x, residual = trial, new_residual
    return x


def canonicalize(mf):
    """Diagonalize the Fock matrix inside the occupied and the virtual blocks separately,
    keeping each orbital's occupation (no aufbau reordering across the blocks)."""
    _, fock, _ = mf.evaluate(mf.density(mf.mo_coeff))
    focks = fock if mf.unrestricted else fock[None]
    coeffs = mf.mo_coeff if mf.unrestricted else mf.mo_coeff[None]
    occs = mf.mo_occ if mf.unrestricted else mf.mo_occ[None]
    new_coeff, new_energy, violation = [], [], []
    for f, c, occ in zip(focks, coeffs, occs):
        c = c.copy()
        energies = np.empty(c.shape[1])
        blocks = [np.flatnonzero(occ > 1e-8), np.flatnonzero(occ <= 1e-8)]
        for block in blocks:
            if len(block):
                e, u = eigh(c[:, block].T @ f @ c[:, block])
                c[:, block] = c[:, block] @ u
                energies[block] = e
        new_coeff.append(c)
        new_energy.append(energies)
        if len(blocks[0]) and len(blocks[1]):
            violation.append(float(energies[blocks[0]].max() - energies[blocks[1]].min()))
    mf.mo_coeff = np.asarray(new_coeff) if mf.unrestricted else new_coeff[0]
    mf.mo_energy = np.asarray(new_energy) if mf.unrestricted else new_energy[0]
    # Positive: the highest occupied orbital lies above the lowest virtual one.
    mf.aufbau_violation_hartree = max(violation) if violation else None


def newton_kernel(mf, dm0=None):
    config = mf.config
    if dm0 is not None:  # Orbitals of the Fock matrix built from the given density.
        _, fock, _ = mf.evaluate(dm0)
        mf.mo_energy, mf.mo_coeff = mf.diagonalize(fock)
    radius, energy, previous_dm = 0.2, None, mf.density(mf.mo_coeff)
    mf.converged = False
    mf.newton_escapes = []
    for cycle in range(mf.max_cycle):
        check_cancelled(mf.directory)
        rotations = OrbitalRotations(mf)
        current, gradient = rotations.energy_and_gradient(np.zeros(rotations.dimension))
        # g = 2 n F_ov; compare the occupied-virtual Fock block with the native residual tolerance.
        fock_ov = np.linalg.norm(gradient) / (2 * rotations.occupation)
        delta = None if energy is None else current - energy
        energy = current
        dm = mf.density(mf.mo_coeff)
        change = float(np.linalg.norm(dm - previous_dm))
        previous_dm = dm
        mf._dm = dm
        mf.callback({'cycle': cycle, 'e_tot': current, 'norm_gorb': fock_ov, 'norm_ddm': change,
                     'delta_energy': delta, 'dm': dm})
        if (fock_ov < config.native_residual_tol / 2 and delta is not None and abs(delta) < config.conv_tol
                and change < config.native_density_tol):
            # A vanishing gradient also holds at saddle points, which Newton steps cannot leave
            # on their own. Accept only a local minimum: otherwise step down the negative mode.
            curvature, mode = lowest_mode(rotations)
            mf.newton_escapes.append({'cycle': cycle + 1, 'lowest_eigenvalue_hartree': curvature})
            if curvature > -1e-4:
                mf.converged = True
                break
            trials = [(rotations.energy_and_gradient(t * mode)[0], t) for t in (.05, .1, .2, .4, .8, -.1, -.4)]
            best, t = min(trials)
            if best >= current:
                break  # No descent along the unstable mode: report as not converged.
            mf.mo_coeff = np.asarray(rotations.rotated(t * mode)) if mf.unrestricted else rotations.rotated(t * mode)[0]
            energy, radius = None, 0.2
            continue
        while True:
            gradient_norm = np.linalg.norm(gradient)
            step = steihaug(rotations.hessian_vector, gradient, radius,
                            tolerance=min(0.1, np.sqrt(gradient_norm)) * gradient_norm,
                            max_iterations=max(20, min(200, rotations.dimension)),
                            diagonal=rotations.diagonal())
            predicted = gradient @ step + 0.5 * step @ rotations.hessian_vector(step)
            trial, trial_gradient = rotations.energy_and_gradient(step)
            ratio = (trial - current) / predicted if predicted < 0 else -1.0
            # Only bypass an unresolved energy ratio at floating-point resolution,
            # and require the gradient to improve. conv_tol is not an acceptance
            # slack: using it here admits cycling with a significant residual.
            rounding = 16 * np.finfo(float).eps * max(1.0, abs(current))
            tiny_improvement = (abs(trial - current) <= rounding and
                                np.linalg.norm(trial_gradient) < gradient_norm)
            if ratio > 0.1 or tiny_improvement:
                mf.mo_coeff = np.asarray(rotations.rotated(step)) if mf.unrestricted else rotations.rotated(step)[0]
                if ratio > 0.75 and np.linalg.norm(step) > 0.9 * radius:
                    radius = min(2 * radius, 1.0)
                elif ratio < 0.25:
                    radius = max(radius / 4, 1e-8)
                break
            radius /= 4
            if radius < 1e-6:
                break
    canonicalize(mf)
    dm = mf.density(mf.mo_coeff)
    energy, fock, components = mf.evaluate(dm)
    mf._dm = dm
    mf.components = components
    mf.final_residual = float(np.linalg.norm(mf.residual(fock, dm)))
    mf.converged = bool(mf.converged and mf.final_residual < config.native_residual_tol
                        and abs(energy - current) < config.conv_tol
                        and np.linalg.norm(dm - previous_dm) < config.native_density_tol)
    mf.scf_summary = {'e1': components['kinetic'] + components['nuclear_attraction'],
                      'coul': components['coulomb'],
                      'exc': components['xc_semilocal'] + components['exact_exchange'],
                      'nuc': components['nuclear_repulsion']}
    mf.e_tot = energy
    return energy
