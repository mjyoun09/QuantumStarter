"""중성 He 원자의 구대칭 Kohn–Sham LDA SCF (방사 격자). QuantumStarter 자체 구현.

원자 단위(에너지 Hartree, 거리 bohr). 고정 점핵 Z=2, 비상대론, 전자 2개, 1s 이중 점유,
스핀 비분극(제한 닫힌 껍질 singlet), 구대칭. 축약 방사함수 u(r) = r R(r),
전자밀도 n(r) = 2|u(r)|^2 / (4 pi r^2), 전자 수 N = ∫ n(r) 4 pi r^2 dr.

이 모듈에서 직접 구현하는 것: 밀도 구성, 하트리 퍼텐셜(구 안쪽·바깥쪽 전하 모두),
LDA 식(Slater 교환 + VWN5 상관 — Libxc LDA_X / LDA_C_VWN과 같은 정의이며
test_atomic_dft에서 Libxc와 대조), SCF 반복·선형 혼합·수렴 판정, 총에너지(이중 계산 보정).
SciPy는 solver.lowest_radial_states의 삼중대각 고유값 계산에만 쓴다(수소형 검증기와 같은
균일 격자·2차 중앙차분). PySCF는 이 모듈에서 쓰지 않는다(scripts/validate_he_lda.py 비교 전용).
"""
import csv
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from numbers import Integral, Real
from pathlib import Path
import platform

import numpy as np
import scipy

from solver import HARTREE_TO_EV, BOHR_TO_ANGSTROM, exact_reduced_radial, lowest_radial_states, radial_grid
from workflow import VERSION

MODEL = 'he-lda'
ELEMENT = 'He'
NUCLEAR_CHARGE = 2
ELECTRONS = 2
OCCUPATION_1S = 2.0
# 이 값 이하의 밀도[bohr^-3]에서는 eps_xc = v_xc = 0으로 둔다(n=0에서 rs가 무한대가 되는 것 방지).
# 기여는 1e-30 * (4 pi r_max^2) * r_max 수준으로 무시할 수 있다.
DENSITY_FLOOR = 1e-30

XC_FUNCTIONALS = {
    'exchange': {'libxc_name': 'LDA_X', 'libxc_id': 1,
                 'description': 'Slater/Dirac exchange, alpha = 2/3, spin-unpolarized',
                 'reference': 'P. A. M. Dirac, Proc. Camb. Phil. Soc. 26, 376 (1930)'},
    'correlation': {'libxc_name': 'LDA_C_VWN', 'libxc_id': 7,
                    'description': 'Vosko-Wilk-Nusair fit 5 (VWN5) of Ceperley-Alder, paramagnetic',
                    'reference': 'S. H. Vosko, L. Wilk, M. Nusair, Can. J. Phys. 58, 1200 (1980)'},
    'pyscf_xc_string': 'LDA_X,LDA_C_VWN',
    'spin': 'unpolarized',
    'evaluation': 'closed-form formulas in atomic_dft.py (no Libxc call at run time); '
                  'values checked against Libxc 7.0.0 in test_atomic_dft',
}
# Libxc lda_c_vwn 상자성 매개변수(A는 Hartree 단위, 원 논문 Rydberg 값의 1/2).
VWN5_PARAMAGNETIC = {'A': 0.0310907, 'b': 3.72744, 'c': 12.9352, 'x0': -0.10498}


@dataclass(frozen=True)
class AtomicDFTConfig:
    """방사 격자와 SCF 설정. 기본 dr = 30/30000 = 0.001 bohr.

    grid_points는 경계 제외 내부 점 수(solver.Config와 같은 정의). density_tol_electrons는
    ∫|n_out - n_in| 4 pi r^2 dr 기준이며, 고유벡터 정밀도 때문에 약 0.4*eps*2/dr^2
    (dr=0.001에서 ~2e-10 e) 아래로는 내려가지 않는다(quality_warnings 참고).
    """
    r_max_bohr: float = 30.0
    grid_points: int = 29999
    mixing: float = 0.5
    max_iterations: int = 100
    energy_tol_hartree: float = 1e-10
    density_tol_electrons: float = 1e-7
    initial_zeta: float = 27.0 / 16.0

    def __post_init__(self):
        for name, (minimum, maximum) in {'grid_points': (100, 1_000_000),
                                         'max_iterations': (1, 10_000)}.items():
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise ValueError(f'{name}: 정수를 입력하세요.')
            if not minimum <= value <= maximum:
                raise ValueError(f'{name}: {minimum}~{maximum} 범위로 입력하세요.')
        for name, maximum in {'r_max_bohr': 1000.0, 'mixing': 1.0, 'energy_tol_hartree': 1.0,
                              'density_tol_electrons': 2.0, 'initial_zeta': 100.0}.items():
            value = getattr(self, name)
            if (isinstance(value, bool) or not isinstance(value, Real)
                    or not np.isfinite(value) or not 0 < value <= maximum):
                raise ValueError(f'{name}: 0보다 크고 {maximum:g} 이하인 유한한 수를 입력하세요.')

    @property
    def dr_bohr(self):
        return self.r_max_bohr / (self.grid_points + 1)


def load_atomic_config(path) -> AtomicDFTConfig:
    data = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(data, dict):
        raise ValueError('He LDA 설정 파일에는 JSON 객체가 필요합니다.')
    if data.get('model', MODEL) != MODEL:
        raise ValueError(f"model이 '{MODEL}'인 설정 파일이 필요합니다.")
    data = {key: value for key, value in data.items() if key != 'model'}
    unknown = set(data) - set(AtomicDFTConfig.__dataclass_fields__)
    if unknown:
        raise ValueError('알 수 없는 He LDA 설정: ' + ', '.join(sorted(unknown)))
    return AtomicDFTConfig(**data)


def save_atomic_config(config: AtomicDFTConfig, path):
    Path(path).write_text(json.dumps({'model': MODEL, **asdict(config)}, ensure_ascii=False, indent=2),
                          encoding='utf-8')


# ---------------------------------------------------------------- LDA 교환·상관

def slater_exchange(n):
    """Slater(Dirac) 교환: (전자당 에너지 eps_x, 퍼텐셜 v_x = d(n eps_x)/dn) [Hartree]."""
    n = np.asarray(n, dtype=float)
    eps = -0.75 * np.cbrt(3.0 / np.pi) * np.cbrt(n)
    return eps, (4.0 / 3.0) * eps


def vwn5_correlation(n):
    """VWN5 상자성 상관: (eps_c, v_c = eps_c - (rs/3) d eps_c/d rs) [Hartree]. x = sqrt(rs)."""
    n = np.asarray(n, dtype=float)
    a, b, c, x0 = (VWN5_PARAMAGNETIC[key] for key in ('A', 'b', 'c', 'x0'))
    x = np.sqrt(np.cbrt(3.0 / (4.0 * np.pi * n)))
    big_x = x * x + b * x + c
    big_x0 = x0 * x0 + b * x0 + c
    q = np.sqrt(4.0 * c - b * b)
    angle = np.arctan(q / (2.0 * x + b))
    eps = a * (np.log(x * x / big_x) + 2.0 * b / q * angle
               - b * x0 / big_x0 * (np.log((x - x0) ** 2 / big_x) + 2.0 * (b + 2.0 * x0) / q * angle))
    # d/dx atan(q/(2x+b)) = -q/(2X) 이므로 아래 식이 eps를 x로 미분한 것이다.
    deps_dx = a * (2.0 / x - 2.0 * (x + b) / big_x
                   - b * x0 / big_x0 * (2.0 / (x - x0) - 2.0 * (x + b + x0) / big_x))
    return eps, eps - x * deps_dx / 6.0


def lda_xc(n):
    """밀도 n[bohr^-3]에서 교환·상관 성분. DENSITY_FLOOR 이하는 0."""
    n = np.asarray(n, dtype=float)
    terms = {key: np.zeros_like(n) for key in ('eps_x', 'v_x', 'eps_c', 'v_c')}
    mask = n > DENSITY_FLOOR
    terms['eps_x'][mask], terms['v_x'][mask] = slater_exchange(n[mask])
    terms['eps_c'][mask], terms['v_c'][mask] = vwn5_correlation(n[mask])
    terms['eps_xc'] = terms['eps_x'] + terms['eps_c']
    terms['v_xc'] = terms['v_x'] + terms['v_c']
    return terms


# ---------------------------------------------------------------- 밀도와 하트리

def density_from_orbital(u_inner, r_inner, occupation=OCCUPATION_1S):
    """내부 점의 u(r)에서 (방사 밀도 4 pi r^2 n [e/bohr], 밀도 n [bohr^-3])."""
    radial = occupation * u_inner ** 2
    return radial, radial / (4.0 * np.pi * r_inner ** 2)


def electron_count(n, r_inner, dr):
    """N = ∫ n(r) 4 pi r^2 dr. u(0)=u(r_max)=0이므로 사다리꼴 = 내부 점 합 * dr."""
    return float(dr * np.sum(4.0 * np.pi * r_inner ** 2 * n))


def hartree_potential(radial_density, r_inner, dr):
    """v_H(r) = (1/r) ∫_0^r rho dr' + ∫_r^r_max rho/r' dr', rho = 4 pi r^2 n.

    첫 항은 반지름 r 안쪽 전하, 둘째 항은 바깥쪽 구껍질의 기여다. rho(0)=rho(r_max)=0인
    사다리꼴 누적합을 쓴다. 이때 이산 핵은 dr/max(r_i, r_j)로 대칭이므로 v_H는
    이산 E_H = (1/2) sum v_H rho dr의 정확한 미분이다.
    """
    inner_charge = dr * (np.cumsum(radial_density) - 0.5 * radial_density)
    shell = radial_density / r_inner
    outer = dr * (np.cumsum(shell[::-1])[::-1] - 0.5 * shell)
    return inner_charge / r_inner + outer


def kinetic_energy(u_inner, dr, occupation=OCCUPATION_1S):
    """T_s = f <u| -1/2 d^2/dr^2 |u>, 해밀토니안과 같은 3점 차분과 u(0)=u(r_max)=0.

    부분합으로 바꾼 같은 값 (f/2) sum_i (u_{i+1}-u_i)^2/dr를 쓴다. 3점 차분 u_{i+1}-2u_i+u_{i-1}의
    자릿수 상쇄가 없어 반올림 오차가 1e-13 Ha 수준이다.
    """
    steps = np.diff(np.pad(u_inner, 1))
    return float(occupation * 0.5 * np.dot(steps, steps) / dr)


def rayleigh_quotient(u_inner, potential, dr):
    """정규화한 u의 이산 해밀토니안 기댓값 <u|H|u>/<u|u>.

    eigh_tridiagonal 고유값의 절대 정밀도는 기본 허용오차 eps*||H||_1 ~ eps*2/dr^2
    (dr=0.001에서 ~4e-10 Ha)로 제한된다. 고유벡터 오차는 Rayleigh 몫에 2차로만 들어가므로
    이 값을 1s 고윳값으로 쓰고, 고유값 해법의 원래 값은 진단용으로 기록한다.
    """
    norm = float(dr * np.dot(u_inner, u_inner))
    return (kinetic_energy(u_inner, dr, 1.0) + float(dr * np.dot(potential, u_inner ** 2))) / norm


def density_functional_terms(radial_density, n, r_inner, dr, v_ext):
    """밀도만의 함수인 에너지 성분 [Hartree]과 그 밀도의 하트리·XC 퍼텐셜."""
    v_hartree = hartree_potential(radial_density, r_inner, dr)
    xc = lda_xc(n)
    terms = {
        'electron_nuclear': float(dr * np.dot(v_ext, radial_density)),
        'hartree': float(0.5 * dr * np.dot(v_hartree, radial_density)),
        'exchange': float(dr * np.dot(xc['eps_x'], radial_density)),
        'correlation': float(dr * np.dot(xc['eps_c'], radial_density)),
    }
    terms['exchange_correlation'] = terms['exchange'] + terms['correlation']
    return terms, v_hartree, xc


# ---------------------------------------------------------------- SCF

@dataclass
class AtomicDFTResult:
    config: AtomicDFTConfig
    dr: float
    r_bohr: np.ndarray            # 내부 점 r_i = i dr, i = 1..N
    u: np.ndarray                 # 마지막 1s 축약 방사함수, sum(u^2) dr = 1, u > 0
    density: np.ndarray           # 마지막 출력 밀도 n_out [bohr^-3]
    potentials: dict              # 마지막 반복의 입력 퍼텐셜 성분 [Hartree]
    eigenvalue_1s: float
    total_energy: float
    components: dict
    history: list
    converged: bool
    status: str
    electron_count: float
    initial_guess: dict = field(default_factory=dict)

    @property
    def radial_density(self):
        return 4.0 * np.pi * self.r_bohr ** 2 * self.density


def initial_orbital(r_inner, dr, zeta):
    """정규화한 수소형 1s u(r) = 2 zeta^(3/2) r exp(-zeta r)를 격자에서 다시 정규화한다."""
    u = exact_reduced_radial(r_inner, 1, 0, zeta)
    norm = float(dr * np.dot(u, u))
    return u / np.sqrt(norm), norm


def run_scf(config: AtomicDFTConfig = AtomicDFTConfig(), progress=None) -> AtomicDFTResult:
    """선형 밀도 혼합 SCF. progress(record)는 반복마다 호출된다.

    반복 k: n_in → v_eff = -Z/r + v_H[n_in] + v_xc[n_in] → 최저 s 상태(ε, u) → n_out.
    잔차는 혼합 전 ∫|n_out - n_in| 4 pi r^2 dr [전자]. 에너지는 Kohn–Sham 함수 E[n_out]을
    고윳값 합으로 쓰고 이중 계산을 보정한 식이다:
        E = 2ε - ∫(v_H[n_in] + v_xc[n_in]) n_out + E_H[n_out] + E_xc[n_out]
    (수렴 시 2ε - E_H - ∫v_xc n + E_xc). |ΔE|와 잔차가 모두 기준 미만이면 수렴이다.
    """
    dr, r = radial_grid(config.r_max_bohr, config.grid_points)
    v_ext = -NUCLEAR_CHARGE / r
    u_in, grid_norm = initial_orbital(r, dr, config.initial_zeta)
    radial_in, n_in = density_from_orbital(u_in, r)
    initial = {'type': 'hydrogenic 1s, normalized on the grid', 'zeta': config.initial_zeta,
               'grid_norm_before_renormalization': grid_norm,
               'electron_count': electron_count(n_in, r, dr)}
    history, previous_energy, converged = [], None, False
    for iteration in range(1, config.max_iterations + 1):
        electrons_in = electron_count(n_in, r, dr)
        v_hartree_in = hartree_potential(radial_in, r, dr)
        xc_in = lda_xc(n_in)
        v_eff = v_ext + v_hartree_in + xc_in['v_xc']
        energies, u_full = lowest_radial_states(v_eff, dr, 1)
        u_out = u_full[1:-1, 0]
        if np.sum(u_out) < 0:  # 전체 부호는 임의이다. 저장·그림용으로만 양수로 맞춘다.
            u_out = -u_out
        eigenvalue = rayleigh_quotient(u_out, v_eff, dr)
        radial_out, n_out = density_from_orbital(u_out, r)
        residual = float(dr * np.sum(np.abs(radial_out - radial_in)))
        terms, _, _ = density_functional_terms(radial_out, n_out, r, dr, v_ext)
        energy = (OCCUPATION_1S * eigenvalue
                  - float(dr * np.dot(v_hartree_in + xc_in['v_xc'], radial_out))
                  + terms['hartree'] + terms['exchange_correlation'])
        delta = None if previous_energy is None else energy - previous_energy
        energy_ok = delta is not None and abs(delta) < config.energy_tol_hartree
        density_ok = residual < config.density_tol_electrons
        record = {'iteration': iteration, 'total_energy_hartree': energy,
                  'delta_energy_hartree': delta, 'density_residual_electrons': residual,
                  'eigenvalue_1s_hartree': eigenvalue,
                  'eigensolver_minus_rayleigh_hartree': float(energies[0]) - eigenvalue,
                  'electrons_in': electrons_in,
                  'energy_criterion_met': energy_ok, 'density_criterion_met': density_ok}
        history.append(record)
        if progress is not None:
            progress(record)
        previous_energy = energy
        if energy_ok and density_ok:
            converged = True
            break
        radial_in = (1.0 - config.mixing) * radial_in + config.mixing * radial_out
        n_in = radial_in / (4.0 * np.pi * r ** 2)

    kinetic = kinetic_energy(u_out, dr)
    direct_total = kinetic + terms['electron_nuclear'] + terms['hartree'] + terms['exchange_correlation']
    components = {
        'kinetic': kinetic, **terms,
        'eigenvalue_sum': OCCUPATION_1S * eigenvalue,
        'double_counting_hartree': -float(dr * np.dot(v_hartree_in, radial_out)) + terms['hartree'],
        'double_counting_xc': -float(dr * np.dot(xc_in['v_xc'], radial_out)) + terms['exchange_correlation'],
        'total_direct': direct_total,
        'total_direct_minus_double_counting': direct_total - energy,
    }
    return AtomicDFTResult(
        config=config, dr=dr, r_bohr=r, u=u_out, density=n_out,
        potentials={'v_ext': v_ext, 'v_hartree': v_hartree_in, 'v_xc': xc_in['v_xc'], 'v_eff': v_eff},
        eigenvalue_1s=eigenvalue, total_energy=energy, components=components, history=history,
        converged=converged, status='converged' if converged else 'not_converged',
        electron_count=electron_count(n_out, r, dr), initial_guess=initial)


# ---------------------------------------------------------------- 결과 저장

def eigensolver_precision(dr):
    """삼중대각 고유값 해법의 기본 절대 허용오차 eps*||H||_1 중 운동 항 eps*2/dr^2 [Hartree]."""
    return float(np.finfo(float).eps * 2.0 / dr ** 2)


def quality_warnings(result: AtomicDFTResult) -> list:
    """검증을 대체하지 않는 간단한 경고."""
    warnings = []
    last = result.history[-1]
    if not result.converged:
        unmet = [name for name, ok in (('|ΔE|', last['energy_criterion_met']),
                                       ('밀도 잔차', last['density_criterion_met'])) if not ok]
        warnings.append(f"SCF 미수렴: 반복 상한 {result.config.max_iterations}회 도달, "
                        f"기준 미달 {', '.join(unmet)}. 마지막 값은 수렴 결과가 아닙니다.")
    if result.config.density_tol_electrons < 5.0 * eigensolver_precision(result.dr):
        warnings.append(f'밀도 잔차 기준이 이 격자의 부동소수점 하한(약 {eigensolver_precision(result.dr):.1e} 수준)에 '
                        '가깝습니다. 기준을 늘리거나 격자를 성기게 하세요.')
    if result.dr * NUCLEAR_CHARGE > 0.01:
        warnings.append('격자가 거칩니다(Z·dr > 0.01). grid_points를 늘려 에너지 수렴을 확인하세요.')
    if result.config.r_max_bohr < 10.0:
        warnings.append('계산영역이 작을 수 있습니다(r_max < 10 bohr). r_max_bohr를 늘려 비교하세요.')
    if abs(result.electron_count - ELECTRONS) > 1e-8:
        warnings.append(f'적분 전자 수가 2와 다릅니다: {result.electron_count:.12f}')
    return warnings


def build_atomic_summary(result: AtomicDFTResult, elapsed_seconds=None) -> dict:
    last = result.history[-1]
    cfg = result.config
    return {
        'program': 'QuantumStarter', 'version': VERSION, 'model': MODEL,
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'description': 'spherical closed-shell Kohn-Sham LDA SCF of the neutral He atom on a radial grid',
        'system': {'element': ELEMENT, 'nuclear_charge': NUCLEAR_CHARGE, 'electrons': ELECTRONS, 'charge': 0,
                   'occupations': {'1s': OCCUPATION_1S}, 'spin': 'unpolarized (restricted singlet)'},
        'assumptions': ['fixed infinite-mass point nucleus', 'nonrelativistic', 'spherical symmetry',
                        'spin-unpolarized LDA', 'Dirichlet radial boundaries u(0)=u(r_max)=0'],
        'method': {
            'grid': 'uniform radial grid r_i = i*dr, dr = r_max/(N+1) (same as solver.py)',
            'kinetic': 'second-order central finite difference, tridiagonal eigensolver (SciPy eigh_tridiagonal)',
            'hartree': 'direct radial integration: inner charge Q(r)/r + outer shells ∫ rho/r dr (trapezoid)',
            'density': 'n(r) = 2|u(r)|^2/(4 pi r^2); electrons = ∫ n 4 pi r^2 dr',
            'mixing': 'linear density mixing n_in <- (1-a) n_in + a n_out',
            'total_energy': 'E = 2*eps_1s - ∫(v_H[n_in]+v_xc[n_in]) n_out + E_H[n_out] + E_xc[n_out]',
        },
        'xc': XC_FUNCTIONALS,
        'config': asdict(cfg),
        'grid': {'dr_bohr': result.dr, 'r_max_bohr': cfg.r_max_bohr, 'interior_points': cfg.grid_points,
                 'r_first_bohr': float(result.r_bohr[0]), 'r_last_bohr': float(result.r_bohr[-1]),
                 'eigensolver_precision_hartree': eigensolver_precision(result.dr)},
        'status': result.status, 'converged': result.converged, 'iterations': len(result.history),
        'convergence': {
            'energy_tol_hartree': cfg.energy_tol_hartree,
            'density_tol_electrons': cfg.density_tol_electrons,
            'final_delta_energy_hartree': last['delta_energy_hartree'],
            'final_density_residual_electrons': last['density_residual_electrons'],
            'energy_criterion_met': last['energy_criterion_met'],
            'density_criterion_met': last['density_criterion_met'],
            'density_residual_definition': '∫|n_out - n_in| 4 pi r^2 dr before mixing [electrons]',
            'delta_energy_definition': 'E(k) - E(k-1) of the total energy below [hartree]',
        },
        'total_energy_hartree': result.total_energy,
        'total_energy_ev': result.total_energy * HARTREE_TO_EV,
        'eigenvalue_1s_hartree': result.eigenvalue_1s,
        'eigenvalue_1s_ev': result.eigenvalue_1s * HARTREE_TO_EV,
        'electron_count': result.electron_count,
        'energy_components_hartree': result.components,
        'initial_guess': result.initial_guess,
        'scf_history': result.history,
        'elapsed_seconds': elapsed_seconds,
        'units': {'energy': 'hartree (eV given separately)', 'distance': 'bohr',
                  'density': 'bohr^-3', 'radial_density': 'electrons/bohr', 'u': 'bohr^(-1/2)'},
        'constants': {'hartree_to_ev': HARTREE_TO_EV, 'bohr_to_angstrom': BOHR_TO_ANGSTROM},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__, 'platform': platform.platform()},
        'warnings': quality_warnings(result),
    }


def save_atomic_result(result: AtomicDFTResult, output_dir, elapsed_seconds=None) -> dict:
    """settings.json, results.json(원자적 교체), scf_history.csv, radial.csv를 저장한다."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = build_atomic_summary(result, elapsed_seconds)
    save_atomic_config(result.config, output_dir / 'settings.json')
    keys = ('iteration', 'total_energy_hartree', 'delta_energy_hartree', 'density_residual_electrons',
            'eigenvalue_1s_hartree', 'electrons_in', 'energy_criterion_met', 'density_criterion_met')
    with (output_dir / 'scf_history.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(keys)
        for record in result.history:
            writer.writerow(['' if record[key] is None else record[key] for key in keys])
    p = result.potentials
    np.savetxt(output_dir / 'radial.csv',
               np.column_stack([result.r_bohr, result.u, result.density, result.radial_density,
                                p['v_ext'], p['v_hartree'], p['v_xc'], p['v_eff']]),
               delimiter=',', comments='',
               header='r_bohr,u_1s,density_bohr^-3,radial_density_4pir2n,v_ext,v_hartree,v_xc,v_eff')
    temporary = output_dir / 'results.json.tmp'
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    temporary.replace(output_dir / 'results.json')
    return summary
