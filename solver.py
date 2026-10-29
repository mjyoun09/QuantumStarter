"""한 전자 원자의 방사형 슈뢰딩거 방정식. 아직 DFT/SCF가 아니다.

원자 단위계: hbar = electron mass = elementary charge = 4*pi*epsilon0 = 1.
u(r) = r*R(r), u(0) = u(r_max) = 0. 고정된 무한 질량 핵을 가정한다.
"""
from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np
from scipy.constants import physical_constants
from scipy.linalg import eigh_tridiagonal
from scipy.special import eval_genlaguerre, gammaln

HARTREE_TO_EV = physical_constants["Hartree energy in eV"][0]
BOHR_TO_ANGSTROM = physical_constants["Bohr radius"][0] / 1e-10


@dataclass(frozen=True)
class Config:
    # 모두 config.json에서 수정할 수 있다. grid_points는 경계 제외 내부 점 수.
    nuclear_charge: int = 1
    angular_momentum_l: int = 0
    num_states: int = 3
    r_max_bohr: float = 60.0
    grid_points: int = 6000

    def __post_init__(self):
        limits = {
            "nuclear_charge": (1, 20),
            "angular_momentum_l": (0, 5),
            "num_states": (1, 6),
            "grid_points": (100, 100000),
        }
        for name, (minimum, maximum) in limits.items():
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise ValueError(f"{name}: 정수를 입력하세요.")
            if not minimum <= value <= maximum:
                raise ValueError(f"{name}: {minimum}~{maximum} 범위로 입력하세요.")
        if (isinstance(self.r_max_bohr, bool)
                or not isinstance(self.r_max_bohr, Real)
                or not np.isfinite(self.r_max_bohr)
                or not 0 < self.r_max_bohr <= 1000):
            raise ValueError("r_max_bohr: 0보다 크고 1000 이하인 유한한 수를 입력하세요.")


@dataclass
class Result:
    config: Config
    r_bohr: np.ndarray          # (N+2,), 두 경계 포함
    u: np.ndarray               # (N+2, states), 축약 방사함수
    u_exact: np.ndarray         # (N+2, states), 무한 공간의 해석해
    energies_hartree: np.ndarray
    exact_energies_hartree: np.ndarray
    principal_n: np.ndarray

    @property
    def dr(self):
        return self.config.r_max_bohr / (self.config.grid_points + 1)

    @property
    def radial_probability(self):
        # P(r) dr = r~r+dr인 구껍질 안에서 전자를 발견할 확률.
        # |psi(x,y,z)|^2인 3차원 부피당 확률밀도와 구별해야 한다.
        return np.abs(self.u) ** 2

    @property
    def mean_r_bohr(self):
        return self.dr * np.sum(self.r_bohr[:, None] * self.radial_probability, axis=0)


def exact_reduced_radial(r_bohr, n, ell, z):
    """정규화된 해석해 u_nl(r). 검증/그래프 전용이며 수치 해에 넣지 않는다."""
    rho = 2.0 * z * np.asarray(r_bohr) / n
    log_norm_sq = (3.0 * np.log(2.0 * z / n)
                   + gammaln(n - ell) - np.log(2.0 * n) - gammaln(n + ell + 1))
    radial_r = (np.exp(0.5 * log_norm_sq) * np.exp(-rho / 2.0)
                * rho ** ell * eval_genlaguerre(n - ell - 1, 2 * ell + 1, rho))
    return np.asarray(r_bohr) * radial_r


def solve(config: Config) -> Result:
    """2차 중앙 유한차분으로 해밀토니안을 만들고 낮은 고유상태를 구한다.

    [-1/2 d^2/dr^2 + l(l+1)/(2r^2) - Z/r] u = E*u
    해석적 에너지 공식은 이 고유값 계산에 사용하지 않는다.
    """
    dr = config.r_max_bohr / (config.grid_points + 1)
    r_inner = dr * np.arange(1, config.grid_points + 1)
    ell = config.angular_momentum_l
    z = config.nuclear_charge

    # r=0은 경계로만 취급하므로 -Z/r의 0 나눗셈이 생기지 않는다.
    potential = -z / r_inner + ell * (ell + 1) / (2.0 * r_inner ** 2)
    diagonal = 1.0 / dr ** 2 + potential
    off_diagonal = np.full(config.grid_points - 1, -0.5 / dr ** 2)

    # SciPy에는 삼중대각 행렬의 고유값 문제만 맡긴다. DFT 엔진은 사용하지 않는다.
    energies, vectors = eigh_tridiagonal(
        diagonal, off_diagonal,
        select="i", select_range=(0, config.num_states - 1),
    )
    # SciPy 정규화 sum(v^2)=1 -> 물리적 정규화 sum(|u|^2)*dr=1.
    u = np.pad(vectors / np.sqrt(dr), ((1, 1), (0, 0)))
    r = dr * np.arange(config.grid_points + 2)
    principal_n = np.arange(ell + 1, ell + config.num_states + 1)
    u_exact = np.column_stack([
        exact_reduced_radial(r, int(n), ell, z) for n in principal_n
    ])
    # 전체 부호는 물리적으로 임의이다. 그림 비교를 위해서만 해석해와 부호를 맞춘다.
    for j in range(config.num_states):
        if np.dot(u[:, j], u_exact[:, j]) < 0:
            u[:, j] *= -1

    return Result(config, r, u, u_exact, energies,
                  -z ** 2 / (2.0 * principal_n.astype(float) ** 2), principal_n)


def quality_warnings(result: Result) -> list[str]:
    """검증을 대체하지 않는 간단한 격자/계산영역 경고."""
    warnings = []
    cfg = result.config
    if result.dr * cfg.nuclear_charge > 0.05:
        warnings.append("격자가 거칩니다. grid_points를 늘려 에너지 수렴을 확인하세요.")
    if cfg.r_max_bohr < 4.0 * int(result.principal_n[-1]) ** 2 / cfg.nuclear_charge:
        warnings.append("계산영역이 작을 수 있습니다. r_max_bohr를 늘려 비교하세요.")
    if np.any(result.energies_hartree >= 0):
        warnings.append("양의 에너지가 나왔습니다. 해당 상태를 원자의 결합상태로 해석하지 마세요.")
    relative_error = np.abs((result.energies_hartree - result.exact_energies_hartree)
                            / result.exact_energies_hartree)
    if np.any(relative_error > 0.01):
        warnings.append("해석해와의 에너지 차이가 1%를 넘습니다. 격자와 경계를 확인하세요.")
    return warnings
