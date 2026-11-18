"""Own PBE (101+130) closed forms and first derivatives; NumPy only at runtime.

PBE: doi:10.1103/PhysRevLett.77.3865, PW92: doi:10.1103/PhysRevB.45.13244.
Constants checked against Libxc 7.0.0 maple/{gga,lda}_exc/{gga_c_pbe,lda_c_pw}.mpl.
PW92 uses the modified (ID 13) coefficients, not ID 12. No Libxc code is called.
All units are atomic. Derivatives refer to energy per volume n*eps, not eps.
"""
import numpy as np

DENSITY_FLOOR = 1e-30
SPIN_FRACTION_FLOOR = 1e-30
BETA = 0.06672455060314922
GAMMA = (1.0 - np.log(2.0)) / np.pi ** 2
MU = BETA * np.pi ** 2 / 3.0
KAPPA = 0.804
PW92 = ((0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294),
        (0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517),
        (0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671))
FZ20 = 1.709920934161365617563962776245
PBE_FUNCTIONALS = {
    'exchange': {'libxc_name': 'GGA_X_PBE', 'libxc_id': 101},
    'correlation': {'libxc_name': 'GGA_C_PBE', 'libxc_id': 130, 'lda_reference': 'LDA_C_PW_MOD (13)'},
    'pyscf_xc_string': 'PBE',
    'implementation': 'own closed forms and first-order chain rule; no runtime Libxc/PySCF',
    'parameters': {'kappa': KAPPA, 'mu': MU, 'beta': BETA, 'gamma': GAMMA},
    'density_floor_bohr_minus3': DENSITY_FLOOR,
    'spin_fraction_floor': SPIN_FRACTION_FLOOR,
    'radial_derivative': 'second-order centered in x; one-sided endpoints; weighted transpose potential',
}


class _First:
    """Small forward chain rule for the 2 or 5 local density/sigma variables.

    Keeps the defining energy and its derivatives together. This is analytic
    differentiation (no finite differences and no external autodiff dependency).
    """
    def __init__(self, value, derivative):
        self.v, self.d = value, derivative

    def __add__(self, other):
        return _First(self.v + other.v, self.d + other.d) if isinstance(other, _First) else _First(self.v + other, self.d)

    __radd__ = __add__

    def __neg__(self):
        return _First(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, _First):
            return _First(self.v * other.v, self.d * other.v + self.v * other.d)
        return _First(self.v * other, self.d * other)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, _First):
            value = self.v / other.v
            return _First(value, (self.d - value * other.d) / other.v)
        return _First(self.v / other, self.d / other)

    def __rtruediv__(self, other):
        return self ** -1 * other

    def __pow__(self, power):
        return _First(self.v ** power, (power * self.v ** (power - 1)) * self.d)

    def log1p(self):
        return _First(np.log1p(self.v), self.d / (1 + self.v))

    def expm1(self):
        return _First(np.expm1(self.v), np.exp(self.v) * self.d)

    def floor(self, minimum):
        return _First(np.maximum(self.v, minimum), self.d * (self.v > minimum))


def _exchange(n, sigma):
    # Zero-spin channels have exactly zero exchange; the internal floor avoids
    # division by zero, and the mask also removes its derivative.
    keep = n.v > DENSITY_FLOOR
    n = n.floor(DENSITY_FLOOR)
    s2 = sigma / (4 * (3 * np.pi ** 2) ** (2 / 3) * n ** (8 / 3))
    enhancement = 1 + KAPPA - KAPPA / (1 + MU / KAPPA * s2)
    return -0.75 * (3 / np.pi) ** (1 / 3) * n ** (4 / 3) * enhancement * keep


def _pw_g(rs, coefficients):
    a, alpha, b1, b2, b3, b4 = coefficients
    q = 2 * a * (b1 * rs ** 0.5 + b2 * rs + b3 * rs ** 1.5 + b4 * rs ** 2)
    return -2 * a * (1 + alpha * rs) * (1 / q).log1p()


def _correlation(n, sigma, up=None, down=None):
    rs = (3 / (4 * np.pi)) ** (1 / 3) * n ** (-1 / 3)
    eps = _pw_g(rs, PW92[0])
    phi = 1.0
    if up is not None:
        plus, minus = 2 * up / n, 2 * down / n
        # The fractional-power derivative at an empty spin is singular. Apply
        # the stated floor to the ENERGY and differentiate that same function.
        # This leaves the occupied-spin limit intact; minority v at exactly zero
        # is a numerical convention, not a finite physical PBE derivative.
        plus, minus = plus.floor(SPIN_FRACTION_FLOOR), minus.floor(SPIN_FRACTION_FLOOR)
        zeta = (up - down) / n
        f = (plus ** (4 / 3) + minus ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
        z4 = zeta ** 4
        eps = eps + (_pw_g(rs, PW92[1]) - eps) * f * z4 - _pw_g(rs, PW92[2]) * f * (1 - z4) / FZ20
        phi = (plus ** (2 / 3) + minus ** (2 / 3)) / 2
    gamma_phi3 = GAMMA * phi ** 3
    y = eps / gamma_phi3
    a = (BETA / GAMMA) / (-y).expm1()
    # t^2 = sigma / [4 phi^2 k_s^2 n^2], k_s^2 = 4 k_F/pi.
    t2 = sigma / ((16 / np.pi) * (3 * np.pi ** 2) ** (1 / 3) * phi ** 2 * n ** (7 / 3))
    u = a * t2
    # Algebraically eps+H, but without subtracting almost equal terms in the
    # large-gradient tail: exp((eps+H)/(gamma phi^3)) = 1+expm1(y)/(1+u+u^2).
    return n * gamma_phi3 * (y.expm1() / (1 + u + u * u)).log1p()


def pbe_pointwise(densities, sigma):
    """Return eps_{x,c,xc}, vrho_* and vsigma_* on independent local variables.

    densities: (1,N) total density or (2,N) spin densities.
    sigma: (1,N) |grad n|^2 or (3,N) [uu,ud,dd]. A physical Gram
    matrix is required. Output derivatives have the corresponding input shapes.
    Below total density 1e-30 all terms are zero. No Libxc 1e-15 clamp is applied.
    """
    densities, sigma = np.asarray(densities, dtype=float), np.asarray(sigma, dtype=float)
    if densities.ndim != 2 or densities.shape[0] not in (1, 2):
        raise ValueError('densities must have shape (1,N) or (2,N)')
    ns, count = densities.shape
    if sigma.shape != ((1 if ns == 1 else 3), count):
        raise ValueError('sigma shape must be (1,N) or (3,N), matching densities')
    if not np.all(np.isfinite(densities)) or np.any(densities < 0) or not np.all(np.isfinite(sigma)):
        raise ValueError('densities must be finite and nonnegative; sigma must be finite')
    if np.any(sigma[0] < 0) or (ns == 2 and (np.any(sigma[2] < 0) or
            np.any(np.abs(sigma[1]) > np.sqrt(np.maximum(sigma[0], np.finfo(float).tiny))
                   * np.sqrt(np.maximum(sigma[2], np.finfo(float).tiny)) * (1 + 1e-12)))):
        raise ValueError('sigma must be a positive semidefinite gradient Gram matrix')
    total = densities.sum(axis=0)
    mask = total > DENSITY_FLOOR
    variables = np.concatenate((densities[:, mask], sigma[:, mask]))
    jets = []
    for index, value in enumerate(variables):
        derivative = np.zeros_like(variables)
        derivative[index] = 1.0
        jets.append(_First(value, derivative))
    if ns == 1:
        n, s = jets
        ex, ec = _exchange(n, s), _correlation(n, s)
    else:
        up, down, suu, sud, sdd = jets
        ex = (_exchange(2 * up, 4 * suu) + _exchange(2 * down, 4 * sdd)) / 2
        ec = _correlation(up + down, suu + 2 * sud + sdd, up, down)
    terms = {}
    for name, energy in (('x', ex), ('c', ec), ('xc', ex + ec)):
        eps, derivative = np.zeros(count), np.zeros((len(jets), count))
        eps[mask], derivative[:, mask] = energy.v / total[mask], energy.d
        terms['eps_' + name] = eps
        terms['vrho_' + name], terms['vsigma_' + name] = derivative[:ns], derivative[ns:]
    return terms


def radial_gradient(density, grid):
    """D n = J^-1 d n/dx; second order including both endpoints."""
    return np.gradient(density, grid.h, axis=-1, edge_order=2) / grid.jac


def radial_gradient_transpose(values, grid):
    """Exact Euclidean transpose of radial_gradient, without a dense matrix."""
    f = values / (2 * grid.h * grid.jac)
    out = np.zeros_like(f)
    out[..., :-2] -= f[..., 1:-1]
    out[..., 2:] += f[..., 1:-1]
    out[..., 0] -= 3 * f[..., 0]
    out[..., 1] += 4 * f[..., 0]
    out[..., 2] -= f[..., 0]
    out[..., -3] += f[..., -1]
    out[..., -2] -= 4 * f[..., -1]
    out[..., -1] += 3 * f[..., -1]
    return out


def pbe_radial(densities, grid):
    """Discrete variational potential v = e_n + W^-1 D^T W e_grad.

    W = 4 pi r^2 J h. The transpose implements -div in spherical coordinates
    with the same boundaries/quadrature as the energy. Do not replace it with
    another call to -radial_gradient: that breaks the discrete chain rule.
    """
    gradient = radial_gradient(densities, grid)
    sigma = gradient ** 2 if len(densities) == 1 else np.array(
        [gradient[0] ** 2, gradient[0] * gradient[1], gradient[1] ** 2])
    terms = pbe_pointwise(densities, sigma)
    vs = terms['vsigma_xc']
    flux = 2 * vs * gradient if len(densities) == 1 else np.array(
        [2 * vs[0] * gradient[0] + vs[1] * gradient[1],
         2 * vs[2] * gradient[1] + vs[1] * gradient[0]])
    weight = 4 * np.pi * grid.r ** 2 * grid.weights
    terms['v_xc'] = terms['vrho_xc'] + radial_gradient_transpose(weight * flux, grid) / weight
    return terms
