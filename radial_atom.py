"""구대칭 원자 Kohn–Sham LDA/LSDA (로그 방사 격자). QuantumStarter 자체 구현.

원자 단위(에너지 Hartree, 거리 bohr). 고정 점핵, 비상대론. 밀도는 구면 평균이다: 부분 껍질은
2l+1개 m 성분을 균등하게 점유한다고 본다(NIST Atomic Reference Data와 같은 관례).
LDA는 스핀 비분극, LSDA는 스핀 분극이며 각 껍질을 다수 스핀부터 2l+1개까지 채운다(Hund 최대 스핀).

격자: r = a (e^x - 1), x_i = i h (i = 1..N). r(0) = 0과 바깥 끝 r(x_{N+1}) ≈ r_max는 Dirichlet 경계.
u(r) = J^(1/2) φ(x), J = dr/dx = r + a로 바꾸면 방사 방정식은 정확히
    -1/2 φ'' + [1/8 + J^2 (l(l+1)/(2 r^2) + v(r))] φ = ε J^2 φ
이다(Liouville 변환). 2차 중앙차분 A φ = ε J^2 φ를 C = J^-1 A J^-1로 대칭화하면 삼중대각이므로
수소형 검증기와 같은 SciPy eigh_tridiagonal을 쓴다. C는 원점 근처 성분이 ~1/(h a)^2인 계단형 행렬이라
기본 허용오차 eps*||C||_1로는 고유값 오차가 Hartree 단위까지 커진다(2026-11-17 실측). 이분법(stebz)에
최소 양수 허용오차를 주면 상대 정밀도로 수렴하며(Rayleigh 몫과 상대 ~1e-12 일치), ε는 상쇄 없는
형태의 Rayleigh 몫으로 다시 계산한다. 적분은 x에 대한 사다리꼴(가중 J h)이다.

직접 구현: 배치 해석, 초기 밀도(Slater 차폐 수소형), 밀도, 하트리(구 안쪽·바깥쪽 전하), LDA/LSDA 식,
SCF·선형 혼합·수렴 판정, 이중 계산을 보정한 총에너지. 비분극 XC는 0.9.5에서 Libxc와 대조한
atomic_dft의 식을 그대로 쓴다. PySCF는 쓰지 않는다(scripts/validate_atoms.py 비교 전용).
"""
import csv
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from numbers import Integral, Real
from pathlib import Path
import platform
import re

import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal

from atomic_dft import DENSITY_FLOOR, VWN5_PARAMAGNETIC, lda_xc, slater_exchange
from solver import BOHR_TO_ANGSTROM, HARTREE_TO_EV, exact_reduced_radial
from workflow import VERSION

MODEL = 'atom'
SYMBOLS = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
           'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
ORBITAL_LETTERS = 'spdf'
NOBLE_CORES = {'He': '1s2', 'Ne': '[He] 2s2 2p6', 'Ar': '[Ne] 3s2 3p6', 'Kr': '[Ar] 3d10 4s2 4p6'}
# 중성 원자 바닥 배치. NIST Atomic Reference Data 표 머리글(Z = 1–36)과 같다.
GROUND_CONFIGURATIONS = {
    1: '1s1', 2: '1s2', 3: '[He] 2s1', 4: '[He] 2s2', 5: '[He] 2s2 2p1', 6: '[He] 2s2 2p2',
    7: '[He] 2s2 2p3', 8: '[He] 2s2 2p4', 9: '[He] 2s2 2p5', 10: '[He] 2s2 2p6', 11: '[Ne] 3s1',
    12: '[Ne] 3s2', 13: '[Ne] 3s2 3p1', 14: '[Ne] 3s2 3p2', 15: '[Ne] 3s2 3p3', 16: '[Ne] 3s2 3p4',
    17: '[Ne] 3s2 3p5', 18: '[Ne] 3s2 3p6', 19: '[Ar] 4s1', 20: '[Ar] 4s2', 21: '[Ar] 3d1 4s2',
    22: '[Ar] 3d2 4s2', 23: '[Ar] 3d3 4s2', 24: '[Ar] 3d5 4s1', 25: '[Ar] 3d5 4s2', 26: '[Ar] 3d6 4s2',
    27: '[Ar] 3d7 4s2', 28: '[Ar] 3d8 4s2', 29: '[Ar] 3d10 4s1', 30: '[Ar] 3d10 4s2',
    31: '[Ar] 3d10 4s2 4p1', 32: '[Ar] 3d10 4s2 4p2', 33: '[Ar] 3d10 4s2 4p3', 34: '[Ar] 3d10 4s2 4p4',
    35: '[Ar] 3d10 4s2 4p5', 36: '[Ar] 3d10 4s2 4p6'}
# eigh_tridiagonal(stebz) 절대 허용오차. 0 이하이면 SciPy가 eps*||C||_1을 쓰므로 최소 양수의 몇 배를 준다.
EIGEN_TOLERANCE = 4.0 * np.finfo(float).tiny
# Libxc lda_c_vwn(ID 7): 상자성·강자성 Ceperley–Alder 맞춤과 스핀 강성 α_c (VWN 1980).
VWN5_FERROMAGNETIC = {'A': 0.01554535, 'b': 7.06042, 'c': 18.0578, 'x0': -0.32500}
VWN5_SPIN_STIFFNESS = {'A': -1.0 / (6.0 * np.pi ** 2), 'b': 1.13107, 'c': 13.0045, 'x0': -0.0047584}
FZETA_DENOMINATOR = 2.0 ** (4.0 / 3.0) - 2.0
FZETA_SECOND_DERIVATIVE_0 = 4.0 / (9.0 * (2.0 ** (1.0 / 3.0) - 1.0))

XC_FUNCTIONALS = {
    'exchange': {'libxc_name': 'LDA_X', 'libxc_id': 1, 'description': 'Slater/Dirac exchange, alpha = 2/3'},
    'correlation': {'libxc_name': 'LDA_C_VWN', 'libxc_id': 7,
                    'description': 'Vosko-Wilk-Nusair fit 5 (VWN5); LSDA uses the VWN spin interpolation '
                                   'with the spin stiffness alpha_c (same as Libxc LDA_C_VWN)',
                    'reference': 'S. H. Vosko, L. Wilk, M. Nusair, Can. J. Phys. 58, 1200 (1980)'},
    'pyscf_xc_string': 'LDA_X,LDA_C_VWN',
    'evaluation': 'closed-form formulas (atomic_dft.py unpolarized, radial_atom.py spin-polarized); '
                  'no Libxc call at run time; checked against Libxc in test_atomic_dft / test_radial_atom',
}


# ---------------------------------------------------------------- 배치

@dataclass(frozen=True)
class Subshell:
    n: int
    l: int
    occupation: float

    @property
    def label(self):
        return f'{self.n}{ORBITAL_LETTERS[self.l]}'

    @property
    def capacity(self):
        return 2 * (2 * self.l + 1)


def element_number(element):
    if isinstance(element, str) and element.strip().capitalize() in SYMBOLS:
        return SYMBOLS.index(element.strip().capitalize()) + 1
    raise ValueError(f'원소 기호를 확인하세요: {element!r} (지원 H–Kr)')


def parse_configuration(text):
    """'[Ne] 3s2 3p6', '1s^2 2s^2 2p^1', '2p0.5' 같은 표기를 Subshell 목록(n, l 순)으로."""
    tokens, subshells = str(text).replace(',', ' ').split(), {}
    while tokens:
        token = tokens.pop(0)
        core = re.fullmatch(r'\[([A-Za-z]+)\]', token)
        if core:
            name = core.group(1).capitalize()
            if name not in NOBLE_CORES:
                raise ValueError(f'지원하지 않는 코어 표기: {token} ([He], [Ne], [Ar], [Kr])')
            tokens = NOBLE_CORES[name].split() + tokens
            continue
        match = re.fullmatch(r'(\d+)([spdfSPDF])\^?(\d+(?:\.\d+)?)?', token)
        if not match:
            raise ValueError(f'배치 표기를 읽을 수 없습니다: {token!r} (예: [Ne] 3s2 3p6)')
        n, l = int(match.group(1)), ORBITAL_LETTERS.index(match.group(2).lower())
        occupation = float(match.group(3)) if match.group(3) is not None else 1.0
        if n <= l:
            raise ValueError(f'{token}: 주양자수 n은 l보다 커야 합니다.')
        if (n, l) in subshells:
            raise ValueError(f'{token}: 같은 껍질이 두 번 지정되었습니다.')
        if not 0.0 <= occupation <= 2 * (2 * l + 1):
            raise ValueError(f'{token}: 점유수는 0~{2 * (2 * l + 1)}입니다.')
        subshells[(n, l)] = Subshell(n, l, occupation)
    if not subshells:
        raise ValueError('배치가 비어 있습니다.')
    return [subshells[key] for key in sorted(subshells)]


def spin_occupations(subshell, spin_polarized):
    """(스핀 이름, 점유수) 목록. LSDA는 다수 스핀에 min(q, 2l+1)개, 나머지를 소수 스핀에."""
    if not spin_polarized:
        return [('unpolarized', subshell.occupation)]
    majority = min(subshell.occupation, 2 * subshell.l + 1.0)
    return [('majority', majority), ('minority', subshell.occupation - majority)]


def slater_effective_charge(z, subshells, target):
    """초기 추정용 Slater 차폐 유효 핵전하."""
    shield = 0.0
    for other in subshells:
        count = other.occupation - (1.0 if other is target else 0.0)
        if count <= 0:
            continue
        if target.l <= 1:
            if other.n == target.n and other.l <= 1:
                shield += (0.30 if target.n == 1 else 0.35) * count
            elif other.n == target.n - 1:
                shield += 0.85 * count
            elif other.n < target.n - 1:
                shield += count
        else:
            if (other.n, other.l) == (target.n, target.l):
                shield += 0.35 * count
            elif other.n < target.n or (other.n == target.n and other.l < target.l):
                shield += count
    return max(z - shield, 0.5)


# ---------------------------------------------------------------- 설정

@dataclass(frozen=True)
class AtomConfig:
    """원자 계산 설정. configuration이 빈 문자열이면 GROUND_CONFIGURATIONS(NIST 바닥 배치)를 쓴다.

    grid_step은 x 간격 h, grid_scale은 a·Z [bohr](원점 근처 간격 ≈ a·h).
    density_tol_electrons는 스핀별 ∫|ρ_out - ρ_in| dr 합 [전자].
    """
    element: str = 'He'
    configuration: str = ''
    spin_polarized: bool = False
    r_max_bohr: float = 50.0
    grid_step: float = 0.001
    grid_scale: float = 1e-4
    mixing: float = 0.3
    max_iterations: int = 300
    energy_tol_hartree: float = 1e-10
    density_tol_electrons: float = 1e-7

    def __post_init__(self):
        z = element_number(self.element)
        object.__setattr__(self, 'element', SYMBOLS[z - 1])
        if not isinstance(self.configuration, str):
            raise ValueError('configuration: 문자열을 입력하세요.')
        if not isinstance(self.spin_polarized, bool):
            raise ValueError('spin_polarized: true 또는 false를 입력하세요.')
        electrons = sum(s.occupation for s in self.subshells)
        if electrons <= 0 or electrons > z + 1e-12:
            raise ValueError(f'전자 수 {electrons:g}: 0보다 크고 핵전하 Z = {z} 이하여야 합니다(음이온은 지원하지 않음).')
        for name, (minimum, maximum) in {'max_iterations': (1, 100_000)}.items():
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or not minimum <= value <= maximum:
                raise ValueError(f'{name}: {minimum}~{maximum} 범위의 정수를 입력하세요.')
        for name, (minimum, maximum) in {'r_max_bohr': (1.0, 1000.0), 'grid_step': (1e-5, 0.05),
                                         'grid_scale': (1e-8, 1e-1), 'mixing': (0.0, 1.0),
                                         'energy_tol_hartree': (0.0, 1.0),
                                         'density_tol_electrons': (0.0, 10.0)}.items():
            value = getattr(self, name)
            valid = (not isinstance(value, bool) and isinstance(value, Real) and np.isfinite(value)
                     and minimum < value <= maximum)
            if not valid:
                raise ValueError(f'{name}: {minimum:g} 초과 {maximum:g} 이하의 유한한 수를 입력하세요.')
        if self.grid_points > 2_000_000:
            raise ValueError('격자점이 너무 많습니다. grid_step을 늘리세요.')

    @property
    def nuclear_charge(self):
        return element_number(self.element)

    @property
    def configuration_text(self):
        return self.configuration.strip() or GROUND_CONFIGURATIONS[self.nuclear_charge]

    @property
    def subshells(self):
        return parse_configuration(self.configuration_text)

    @property
    def scale_bohr(self):
        return self.grid_scale / self.nuclear_charge

    @property
    def grid_points(self):
        return max(int(round(np.log(self.r_max_bohr / self.scale_bohr + 1.0) / self.grid_step)) - 1, 10)


def load_atom_config(path) -> AtomConfig:
    data = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(data, dict) or data.get('model', MODEL) != MODEL:
        raise ValueError(f"model이 '{MODEL}'인 설정 JSON 객체가 필요합니다.")
    data = {key: value for key, value in data.items() if key != 'model'}
    unknown = set(data) - set(AtomConfig.__dataclass_fields__)
    if unknown:
        raise ValueError('알 수 없는 원자 설정: ' + ', '.join(sorted(unknown)))
    return AtomConfig(**data)


def save_atom_config(config: AtomConfig, path):
    Path(path).write_text(json.dumps({'model': MODEL, **asdict(config)}, ensure_ascii=False, indent=2),
                          encoding='utf-8')


# ---------------------------------------------------------------- 격자와 고유값

@dataclass(frozen=True)
class LogGrid:
    h: float
    a: float
    x: np.ndarray       # 내부 점 x_i = i h
    r: np.ndarray       # r_i = a (e^x_i - 1)
    jac: np.ndarray     # J_i = dr/dx = r_i + a
    r_end: float        # 바깥 Dirichlet 경계 r(x_{N+1})

    @property
    def weights(self):
        """∫ f(r) dr ≈ Σ f_i J_i h (x에 대한 사다리꼴; 양 끝 함수값 0)."""
        return self.jac * self.h


def log_grid(config: AtomConfig) -> LogGrid:
    h, a, count = config.grid_step, config.scale_bohr, config.grid_points
    x = h * np.arange(1, count + 1)
    r = a * np.expm1(x)
    return LogGrid(h, a, x, r, r + a, float(a * np.expm1((count + 1) * h)))


def rayleigh_quotients(grid: LogGrid, phi, transformed_potential):
    """열마다 (ε, 운동에너지) — 정규화 Σ J²φ² h = 1 기준. 차분 제곱 합이라 자릿수 상쇄가 없다."""
    steps = np.diff(np.pad(phi, ((1, 1), (0, 0))), axis=0)
    gradient = 0.5 * np.sum(steps ** 2, axis=0) / grid.h
    norm = grid.h * np.sum((grid.jac[:, None] * phi) ** 2, axis=0)
    potential = grid.h * np.sum(transformed_potential[:, None] * phi ** 2, axis=0)
    return (gradient + potential) / norm, norm


def solve_channel(grid: LogGrid, potential, ell, count):
    """각운동량 ell 채널의 낮은 count개 상태. 반환: (ε Rayleigh, ε 고유값 해법, 운동에너지, φ, 마디 수)."""
    centrifugal = ell * (ell + 1) / (2.0 * grid.r ** 2)
    kinetic_shift = 0.125 + grid.jac ** 2 * centrifugal
    transformed = kinetic_shift + grid.jac ** 2 * potential
    diagonal = (1.0 / grid.h ** 2 + transformed) / grid.jac ** 2
    off_diagonal = -0.5 / (grid.h ** 2 * grid.jac[:-1] * grid.jac[1:])
    raw, vectors = eigh_tridiagonal(diagonal, off_diagonal, select='i', select_range=(0, count - 1),
                                    lapack_driver='stebz', tol=EIGEN_TOLERANCE)
    phi = vectors / (np.sqrt(grid.h) * grid.jac[:, None])            # Σ J²φ² h = 1
    for k in range(count):                                            # 원점 쪽 첫 고리를 양수로
        column = phi[:, k]
        first = np.argmax(np.abs(column) > 1e-8 * np.max(np.abs(column)))
        if column[first] < 0:
            phi[:, k] = -column
    eigenvalues, _ = rayleigh_quotients(grid, phi, transformed)
    kinetic, _ = rayleigh_quotients(grid, phi, kinetic_shift)
    nodes = []
    for k in range(count):
        column = phi[:, k]
        significant = column[np.abs(column) > 1e-6 * np.max(np.abs(column))]
        nodes.append(int(np.count_nonzero(np.diff(np.sign(significant)) != 0)))
    return eigenvalues, raw, kinetic, phi, nodes


def hartree_potential(radial_density, grid: LogGrid):
    """v_H(r) = (1/r)∫_0^r ρ dr' + ∫_r ρ/r' dr' (ρ = 4πr²n). x에 대한 사다리꼴 누적합.

    이산 핵은 J_j h / max(r_i, r_j)로 대칭이므로 v_H는 이산 E_H = ½ Σ v_H ρ J h의 정확한 미분이다.
    """
    charge = radial_density * grid.weights
    inner = np.cumsum(charge) - 0.5 * charge
    shell = charge / grid.r
    outer = np.cumsum(shell[::-1])[::-1] - 0.5 * shell
    return inner / grid.r + outer


# ---------------------------------------------------------------- LSDA 교환·상관

def _vwn_term(x, parameters):
    """VWN 맞춤 함수 G(x) = A[ln(x²/X) + ...]와 dG/dx (x = sqrt(rs))."""
    a, b, c, x0 = (parameters[key] for key in ('A', 'b', 'c', 'x0'))
    big_x = x * x + b * x + c
    big_x0 = x0 * x0 + b * x0 + c
    q = np.sqrt(4.0 * c - b * b)
    angle = np.arctan(q / (2.0 * x + b))
    value = a * (np.log(x * x / big_x) + 2.0 * b / q * angle
                 - b * x0 / big_x0 * (np.log((x - x0) ** 2 / big_x) + 2.0 * (b + 2.0 * x0) / q * angle))
    slope = a * (2.0 / x - 2.0 * (x + b) / big_x
                 - b * x0 / big_x0 * (2.0 / (x - x0) - 2.0 * (x + b + x0) / big_x))
    return value, slope


def vwn5_spin_correlation(n_up, n_down):
    """VWN5 스핀 분극 상관 (Libxc LDA_C_VWN, spin=polarized와 같은 식).

    ε_c = ε_P + α_c f(ζ)(1-ζ⁴)/f''(0) + (ε_F - ε_P) f(ζ) ζ⁴.
    반환: (ε_c, v_c↑, v_c↓) [Hartree]. 밀도는 모두 DENSITY_FLOOR보다 커야 한다.
    """
    n = n_up + n_down
    zeta = np.clip((n_up - n_down) / n, -1.0, 1.0)
    x = np.sqrt(np.cbrt(3.0 / (4.0 * np.pi * n)))
    para, dpara = _vwn_term(x, VWN5_PARAMAGNETIC)
    ferro, dferro = _vwn_term(x, VWN5_FERROMAGNETIC)
    stiff, dstiff = _vwn_term(x, VWN5_SPIN_STIFFNESS)
    plus, minus = 1.0 + zeta, 1.0 - zeta
    f = (plus ** (4.0 / 3.0) + minus ** (4.0 / 3.0) - 2.0) / FZETA_DENOMINATOR
    df = (4.0 / 3.0) * (np.cbrt(plus) - np.cbrt(minus)) / FZETA_DENOMINATOR
    z3, z4 = zeta ** 3, zeta ** 4
    eps = para + stiff * f * (1.0 - z4) / FZETA_SECOND_DERIVATIVE_0 + (ferro - para) * f * z4
    deps_dx = dpara + dstiff * f * (1.0 - z4) / FZETA_SECOND_DERIVATIVE_0 + (dferro - dpara) * f * z4
    deps_dzeta = (stiff * (df * (1.0 - z4) - 4.0 * z3 * f) / FZETA_SECOND_DERIVATIVE_0
                  + (ferro - para) * (df * z4 + 4.0 * z3 * f))
    common = eps - x * deps_dx / 6.0                                  # ε - (rs/3) ∂ε/∂rs
    return eps, common + minus * deps_dzeta, common - plus * deps_dzeta


def lsda_xc(n_up, n_down):
    """스핀 밀도에서 교환·상관 성분. eps_*는 전체 전자당, v_*는 (2, N) 스핀별 퍼텐셜."""
    n_up, n_down = np.asarray(n_up, dtype=float), np.asarray(n_down, dtype=float)
    n = n_up + n_down
    terms = {'eps_x': np.zeros_like(n), 'eps_c': np.zeros_like(n),
             'v_x': np.zeros((2,) + n.shape), 'v_c': np.zeros((2,) + n.shape)}
    mask = n > DENSITY_FLOOR
    up, down = n_up[mask], n_down[mask]
    ex_up, vx_up = slater_exchange(2.0 * up)
    ex_down, vx_down = slater_exchange(2.0 * down)
    terms['eps_x'][mask] = (up * ex_up + down * ex_down) / (up + down)
    terms['v_x'][0][mask], terms['v_x'][1][mask] = vx_up, vx_down
    terms['eps_c'][mask], terms['v_c'][0][mask], terms['v_c'][1][mask] = vwn5_spin_correlation(up, down)
    terms['eps_xc'] = terms['eps_x'] + terms['eps_c']
    terms['v_xc'] = terms['v_x'] + terms['v_c']
    return terms


def xc_terms(spin_radial, grid: LogGrid):
    """스핀별 방사 밀도 (nspin, N) → XC 성분. 비분극은 atomic_dft.lda_xc를 그대로 쓴다."""
    densities = spin_radial / (4.0 * np.pi * grid.r ** 2)
    if densities.shape[0] == 1:
        terms = lda_xc(densities[0])
        return {**terms, 'v_xc': terms['v_xc'][None, :]}
    return lsda_xc(densities[0], densities[1])


# ---------------------------------------------------------------- SCF

@dataclass
class AtomResult:
    config: AtomConfig
    grid: LogGrid
    spin_labels: tuple
    orbitals: list               # dict: label, n, l, spin, occupation, eigenvalue_hartree, ...
    wavefunctions: dict          # (label, spin) -> u(r) 내부 점, ∫u² dr = 1
    spin_radial: np.ndarray      # (nspin, N) 출력 방사 밀도 4πr²n_σ [전자/bohr]
    potentials: dict             # 마지막 입력 퍼텐셜: v_ext, v_hartree, v_xc (nspin, N)
    total_energy: float
    components: dict
    history: list
    converged: bool
    status: str
    electron_count: float
    initial_guess: dict = field(default_factory=dict)


def _channels(subshells, spin_labels, spin_polarized):
    """(spin 순번, l) -> 필요한 상태 수(최대 n - l)."""
    channels = {}
    for subshell in subshells:
        for spin_index, _ in enumerate(spin_labels):
            key = (spin_index, subshell.l)
            channels[key] = max(channels.get(key, 0), subshell.n - subshell.l)
    return channels


def initial_density(config: AtomConfig, grid: LogGrid, spin_labels):
    z, subshells = config.nuclear_charge, config.subshells
    radial = np.zeros((len(spin_labels), grid.r.size))
    charges = {}
    for subshell in subshells:
        z_eff = slater_effective_charge(z, subshells, subshell)
        charges[subshell.label] = z_eff
        u = exact_reduced_radial(grid.r, subshell.n, subshell.l, z_eff)
        u = u / np.sqrt(np.sum(u * u * grid.weights))
        for spin_index, (_, occupation) in enumerate(spin_occupations(subshell, config.spin_polarized)):
            radial[spin_index] += occupation * u * u
    return radial, {'type': 'hydrogenic orbitals with Slater-rule effective charges, normalized on the grid',
                    'effective_charges': charges}


def run_atom(config: AtomConfig = AtomConfig(), progress=None) -> AtomResult:
    """선형 스핀 밀도 혼합 SCF. progress(record)는 반복마다 호출된다.

    반복 k: ρ_in,σ → v_σ = -Z/r + v_H[ρ_in] + v_xc,σ[ρ_in] → 채널(σ, l)별 최저 상태 → ρ_out,σ.
    잔차는 혼합 전 Σ_σ ∫|ρ_out,σ - ρ_in,σ| dr. 총에너지는 고윳값 합의 이중 계산 보정:
        E = Σ f ε - Σ_σ ∫(v_H + v_xc,σ)[ρ_in] ρ_out,σ dr + E_H[ρ_out] + E_xc[ρ_out]
    |ΔE|와 잔차가 모두 기준 미만이면 수렴이다.
    """
    z, subshells = config.nuclear_charge, config.subshells
    grid = log_grid(config)
    spin_labels = ('majority', 'minority') if config.spin_polarized else ('unpolarized',)
    channels = _channels(subshells, spin_labels, config.spin_polarized)
    v_ext = -z / grid.r
    radial_in, initial = initial_density(config, grid, spin_labels)
    history, previous_energy, converged = [], None, False
    for iteration in range(1, config.max_iterations + 1):
        total_in = radial_in.sum(axis=0)
        v_hartree = hartree_potential(total_in, grid)
        xc_in = xc_terms(radial_in, grid)
        potentials = v_ext + v_hartree + xc_in['v_xc']
        solutions = {(spin, ell): solve_channel(grid, potentials[spin], ell, count)
                     for (spin, ell), count in channels.items()}
        radial_out = np.zeros_like(radial_in)
        orbitals, wavefunctions, band = [], {}, 0.0
        for subshell in subshells:
            for spin, (spin_name, occupation) in enumerate(spin_occupations(subshell, config.spin_polarized)):
                eigenvalues, raw, kinetic, phi, nodes = solutions[(spin, subshell.l)]
                k = subshell.n - subshell.l - 1
                u = np.sqrt(grid.jac) * phi[:, k]
                radial_out[spin] += occupation * u * u
                band += occupation * eigenvalues[k]
                wavefunctions[(subshell.label, spin_name)] = u
                orbitals.append({'label': subshell.label, 'n': subshell.n, 'l': subshell.l, 'spin': spin_name,
                                 'occupation': occupation, 'eigenvalue_hartree': float(eigenvalues[k]),
                                 'eigensolver_minus_rayleigh_hartree': float(raw[k] - eigenvalues[k]),
                                 'kinetic_hartree': float(kinetic[k]), 'nodes': nodes[k],
                                 'expected_nodes': subshell.n - subshell.l - 1})
        residual = float(np.sum(np.abs(radial_out - radial_in) * grid.weights))
        total_out = radial_out.sum(axis=0)
        xc_out = xc_terms(radial_out, grid)
        hartree_out = 0.5 * float(np.sum(hartree_potential(total_out, grid) * total_out * grid.weights))
        exchange = float(np.sum(xc_out['eps_x'] * total_out * grid.weights))
        correlation = float(np.sum(xc_out['eps_c'] * total_out * grid.weights))
        double_counting = float(np.sum((v_hartree + xc_in['v_xc']) * radial_out * grid.weights))
        band = float(band)
        energy = band - double_counting + hartree_out + exchange + correlation
        delta = None if previous_energy is None else energy - previous_energy
        energy_ok = bool(delta is not None and abs(delta) < config.energy_tol_hartree)
        density_ok = bool(residual < config.density_tol_electrons)
        record = {'iteration': iteration, 'total_energy_hartree': energy, 'delta_energy_hartree': delta,
                  'density_residual_electrons': residual, 'energy_criterion_met': energy_ok,
                  'density_criterion_met': density_ok,
                  'highest_occupied_eigenvalue_hartree': max(o['eigenvalue_hartree'] for o in orbitals
                                                             if o['occupation'] > 0)}
        history.append(record)
        if progress is not None:
            progress(record)
        previous_energy = energy
        if energy_ok and density_ok:
            converged = True
            break
        radial_in = (1.0 - config.mixing) * radial_in + config.mixing * radial_out

    kinetic_total = sum(o['occupation'] * o['kinetic_hartree'] for o in orbitals)
    electron_nuclear = float(np.sum(v_ext * total_out * grid.weights))
    direct = kinetic_total + electron_nuclear + hartree_out + exchange + correlation
    components = {'kinetic': kinetic_total, 'electron_nuclear': electron_nuclear, 'hartree': hartree_out,
                  'exchange': exchange, 'correlation': correlation, 'exchange_correlation': exchange + correlation,
                  'eigenvalue_sum': band, 'total_direct': direct, 'total_direct_minus_double_counting': direct - energy}
    return AtomResult(config=config, grid=grid, spin_labels=spin_labels, orbitals=orbitals,
                      wavefunctions=wavefunctions, spin_radial=radial_out,
                      potentials={'v_ext': v_ext, 'v_hartree': v_hartree, 'v_xc': xc_in['v_xc']},
                      total_energy=energy, components=components, history=history, converged=converged,
                      status='converged' if converged else 'not_converged',
                      electron_count=float(np.sum(total_out * grid.weights)), initial_guess=initial)


# ---------------------------------------------------------------- 결과 저장

def quality_warnings(result: AtomResult) -> list:
    warnings, last, cfg = [], result.history[-1], result.config
    if not result.converged:
        unmet = [name for name, ok in (('|ΔE|', last['energy_criterion_met']),
                                       ('밀도 잔차', last['density_criterion_met'])) if not ok]
        warnings.append(f"SCF 미수렴: 반복 상한 {cfg.max_iterations}회 도달, 기준 미달 {', '.join(unmet)}. "
                        '마지막 값은 수렴 결과가 아닙니다.')
    wrong = [f"{o['label']}({o['spin']})" for o in result.orbitals if o['nodes'] != o['expected_nodes']]
    if wrong:
        warnings.append('마디 수가 n-l-1과 다른 오비탈: ' + ', '.join(wrong) + '. 격자를 확인하세요.')
    unbound = [f"{o['label']}({o['spin']})" for o in result.orbitals
               if o['occupation'] > 0 and o['eigenvalue_hartree'] >= 0]
    if unbound:
        warnings.append('점유된 오비탈의 고윳값이 0 이상입니다(속박 안 됨): ' + ', '.join(unbound))
    if cfg.grid_step > 0.005:
        warnings.append('격자가 거칩니다(grid_step > 0.005). 상대 이산화 오차 ~0.05·h²를 확인하세요.')
    electrons = sum(s.occupation for s in cfg.subshells)
    if abs(result.electron_count - electrons) > 1e-8 * max(electrons, 1):
        warnings.append(f'적분 전자 수 {result.electron_count:.12f}가 배치의 {electrons:g}와 다릅니다.')
    return warnings


def build_atom_summary(result: AtomResult, elapsed_seconds=None) -> dict:
    cfg, grid, last = result.config, result.grid, result.history[-1]
    electrons = sum(s.occupation for s in cfg.subshells)
    return {
        'program': 'QuantumStarter', 'version': VERSION, 'model': MODEL,
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'description': 'spherical Kohn-Sham ' + ('LSDA' if cfg.spin_polarized else 'LDA')
                       + ' SCF of an atom on a logarithmic radial grid',
        'system': {'element': cfg.element, 'nuclear_charge': cfg.nuclear_charge, 'electrons': electrons,
                   'charge': cfg.nuclear_charge - electrons, 'configuration': cfg.configuration_text,
                   'spin': ('spin-polarized, each subshell filled majority-first (Hund)' if cfg.spin_polarized
                            else 'unpolarized'),
                   'spin_moment_electrons': (float(sum(o['occupation'] * (1 if o['spin'] == 'majority' else -1)
                                                       for o in result.orbitals)) if cfg.spin_polarized else 0.0)},
        'assumptions': ['fixed infinite-mass point nucleus', 'nonrelativistic',
                        'spherically averaged density (partially filled subshells averaged over m)',
                        'Dirichlet radial boundaries u(0)=u(r_end)=0'],
        'method': {
            'grid': 'r = a(exp(x)-1), uniform x, u = J^(1/2) phi, J = r + a (exact Liouville transform)',
            'eigensolver': 'symmetric tridiagonal C = J^-1 A J^-1, SciPy eigh_tridiagonal stebz with '
                           'tol = 4*tiny; eigenvalue = Rayleigh quotient',
            'hartree': 'inner charge Q(r)/r + outer shells ∫ rho/r dr (trapezoid in x)',
            'mixing': 'linear mixing of spin radial densities',
            'total_energy': 'E = sum f eps - sum_s ∫(v_H + v_xc,s)[in] rho_out,s + E_H[out] + E_xc[out]'},
        'xc': {**XC_FUNCTIONALS, 'spin': 'polarized' if cfg.spin_polarized else 'unpolarized'},
        'config': asdict(cfg),
        'grid': {'grid_step_h': grid.h, 'scale_a_bohr': grid.a, 'interior_points': int(grid.r.size),
                 'r_first_bohr': float(grid.r[0]), 'r_last_bohr': float(grid.r[-1]), 'r_end_bohr': grid.r_end,
                 'first_spacing_bohr': float(grid.r[0])},
        'status': result.status, 'converged': result.converged, 'iterations': len(result.history),
        'convergence': {
            'energy_tol_hartree': cfg.energy_tol_hartree, 'density_tol_electrons': cfg.density_tol_electrons,
            'final_delta_energy_hartree': last['delta_energy_hartree'],
            'final_density_residual_electrons': last['density_residual_electrons'],
            'energy_criterion_met': last['energy_criterion_met'],
            'density_criterion_met': last['density_criterion_met'],
            'density_residual_definition': 'sum over spins of ∫|rho_out - rho_in| dr before mixing [electrons]'},
        'total_energy_hartree': result.total_energy, 'total_energy_ev': result.total_energy * HARTREE_TO_EV,
        'electron_count': result.electron_count,
        'orbitals': result.orbitals,
        'energy_components_hartree': result.components,
        'initial_guess': result.initial_guess,
        'scf_history': result.history,
        'elapsed_seconds': elapsed_seconds,
        'units': {'energy': 'hartree (eV given for the total)', 'distance': 'bohr',
                  'radial_density': 'electrons/bohr', 'u': 'bohr^(-1/2)'},
        'constants': {'hartree_to_ev': HARTREE_TO_EV, 'bohr_to_angstrom': BOHR_TO_ANGSTROM},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__, 'platform': platform.platform()},
        'warnings': quality_warnings(result),
    }


def save_atom_result(result: AtomResult, output_dir, elapsed_seconds=None) -> dict:
    """settings.json, results.json(원자적 교체), scf_history.csv, orbitals.csv, radial.csv."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = build_atom_summary(result, elapsed_seconds)
    save_atom_config(result.config, output_dir / 'settings.json')
    for name, rows, keys in (
            ('scf_history.csv', result.history,
             ('iteration', 'total_energy_hartree', 'delta_energy_hartree', 'density_residual_electrons',
              'highest_occupied_eigenvalue_hartree', 'energy_criterion_met', 'density_criterion_met')),
            ('orbitals.csv', result.orbitals,
             ('label', 'n', 'l', 'spin', 'occupation', 'eigenvalue_hartree', 'kinetic_hartree', 'nodes',
              'eigensolver_minus_rayleigh_hartree'))):
        with (output_dir / name).open('w', newline='', encoding='utf-8-sig') as stream:
            writer = csv.writer(stream)
            writer.writerow(keys)
            for row in rows:
                writer.writerow(['' if row[key] is None else row[key] for key in keys])
    columns = [result.grid.r] + list(result.spin_radial) + list(result.potentials['v_xc'])
    headers = ['r_bohr'] + [f'radial_density_{s}' for s in result.spin_labels] + [f'v_xc_{s}' for s in result.spin_labels]
    columns += [result.potentials['v_ext'], result.potentials['v_hartree']]
    headers += ['v_ext', 'v_hartree']
    for (label, spin), u in result.wavefunctions.items():
        columns.append(u)
        headers.append(f'u_{label}_{spin}')
    np.savetxt(output_dir / 'radial.csv', np.column_stack(columns), delimiter=',', comments='',
               header=','.join(headers))
    temporary = output_dir / 'results.json.tmp'
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    temporary.replace(output_dir / 'results.json')
    return summary
