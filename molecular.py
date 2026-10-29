"""Engine-independent molecular input. Coordinates are never inferred from a name."""
from dataclasses import asdict, dataclass, field
import json
import math
from numbers import Integral, Real
from pathlib import Path
from element_catalog import MOLECULAR_CAPABILITIES, molecular_core_counts

ELEMENTS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe '
            'Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn '
            'Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W '
            'Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf '
            'Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()
ATOMIC_NUMBERS = {name: i + 1 for i, name in enumerate(ELEMENTS)}
# Recognize all elements. Calculation support is checked against actual basis data.
MAX_ATOMIC_NUMBER = 118
SUPPORTED_ELEMENTS = tuple(ELEMENTS[:MAX_ATOMIC_NUMBER])
PERIOD4_ELEMENTS = tuple(ELEMENTS[18:36])
TRANSITION_METALS = frozenset(ELEMENTS[20:30] + ELEMENTS[38:48] + ELEMENTS[56:80] + ELEMENTS[88:112])
DEF2_ECP_ELEMENTS = frozenset(molecular_core_counts('def2-svp'))
DEF2_BASES = frozenset({'def2-svp', 'def2-tzvp'})
DEF2_CORE_ELECTRONS = molecular_core_counts('def2-svp')
BASES = tuple(MOLECULAR_CAPABILITIES)
RELATIVISTIC_BASES = frozenset({'dyall-v2z', 'dyall-v3z'})
FUNCTIONALS = ('lda,vwn', 'pbe', 'b3lyp', 'pbe0')
DISPERSIONS = ('none', 'd4')
CALCULATIONS = ('single_point', 'gradient', 'optimize', 'frequencies', 'scan', 'transition_state', 'neb')
# Thresholds handed to geomeTRIC: energy change [Ha], RMS/max gradient [Ha/bohr], RMS/max step [Angstrom].
# 'default' equals geomeTRIC's own defaults (as documented by PySCF's geomeTRIC driver).
OPT_CONVERGENCE = {
    'default': {'energy': 1e-6, 'grms': 3e-4, 'gmax': 4.5e-4, 'drms': 1.2e-3, 'dmax': 1.8e-3},
    'tight': {'energy': 1e-7, 'grms': 3e-5, 'gmax': 4.5e-5, 'drms': 1.2e-4, 'dmax': 1.8e-4},
}
# Checked against dftd4 4.3.1 (test_dispersion): published D4 damping parameters exist for
# these functionals (none for LDA), and the reference data cover Z=1-103 and 112-118.
D4_FUNCTIONALS = frozenset({'pbe', 'pbe0', 'b3lyp'})
D4_ELEMENTS = frozenset(ELEMENTS[:103] + ELEMENTS[111:118])
# Verified against bundled PySCF basis data. Never silently substitute a basis.
BASIS_UNAVAILABLE = {basis: frozenset(ELEMENTS) - set(info['elements'])
                     for basis, info in MOLECULAR_CAPABILITIES.items()}


def available_bases(atoms):
    elements = {atom.element for atom in atoms}
    return tuple(basis for basis in BASES if not elements.intersection(BASIS_UNAVAILABLE.get(basis, ())))


def parse_basis_overrides(text):
    """'O=def2-tzvp, H=sto-3g' -> {'O': 'def2-tzvp', 'H': 'sto-3g'}; validated later by DFTConfig."""
    overrides = {}
    for item in text.replace(';', ',').replace('\n', ',').split(','):
        item = item.strip()
        if not item:
            continue
        if item.count('=') != 1:
            raise ValueError(f'원소별 기저 "{item}": 원소=기저 형식으로 입력하세요. 예: O=def2-tzvp')
        element, basis = (part.strip() for part in item.split('='))
        element = element.capitalize()
        if element in overrides:
            raise ValueError(f'원소별 기저에 {element}가 두 번 지정되었습니다.')
        overrides[element] = basis.lower()
    return overrides


def format_basis_overrides(overrides):
    return ', '.join(f'{element}={basis}' for element, basis in overrides.items())


def describe_bases(default, basis_by_element):
    """'def2-svp (O=def2-tzvp)': the default plus the elements that differ from it."""
    different = {e: b for e, b in basis_by_element.items() if b != default}
    return default + (f' ({format_basis_overrides(different)})' if different else '')


SCAN_ATOMS = {'distance': 2, 'angle': 3, 'dihedral': 4}


def parse_scan(spec, natoms):
    """'distance 1 2 0.7 1.2 6' -> ('distance', (1, 2), [0.7, ..., 1.2]); Angstrom / degrees, 1-based atoms."""
    parts = spec.split()
    if not parts or parts[0] not in SCAN_ATOMS or len(parts) != SCAN_ATOMS[parts[0]] + 4:
        raise ValueError('스캔 형식: distance i j 시작 끝 점수 / angle i j k ... / dihedral i j k l ... '
                         '(Å·도, 원자 번호는 1부터).')
    count = SCAN_ATOMS[parts[0]]
    try:
        atoms = tuple(int(x) for x in parts[1:1 + count])
        start, stop, points = float(parts[1 + count]), float(parts[2 + count]), int(parts[3 + count])
    except ValueError:
        raise ValueError('스캔의 원자 번호·시작·끝·점수가 숫자가 아닙니다.') from None
    if len(set(atoms)) != count or not all(1 <= a <= natoms for a in atoms):
        raise ValueError(f'스캔 원자 번호는 1~{natoms} 사이의 서로 다른 값이어야 합니다.')
    if not 2 <= points <= 101 or not (math.isfinite(start) and math.isfinite(stop)) or start == stop:
        raise ValueError('스캔 점수는 2~101, 시작과 끝은 서로 다른 유한한 값이어야 합니다.')
    return parts[0], atoms, [start + (stop - start) * k / (points - 1) for k in range(points)]


def integer(value, label, low, high):
    if isinstance(value, bool) or not isinstance(value, Integral) or not low <= value <= high:
        raise ValueError(f'{label}: {low}~{high} 범위의 정수를 입력하세요.')


def real(value, label, low, high):
    if (isinstance(value, bool) or not isinstance(value, Real)
            or not math.isfinite(value) or not low <= value <= high):
        raise ValueError(f'{label}: {low}~{high} 범위의 유한한 수를 입력하세요.')


@dataclass(frozen=True)
class Atom:
    element: str
    position: tuple[float, float, float]

    def __post_init__(self):
        if not isinstance(self.element, str):
            raise ValueError('원소는 H, C, O처럼 기호로 입력하세요.')
        element = self.element.strip().capitalize()
        if element not in ATOMIC_NUMBERS:
            raise ValueError(f'알 수 없는 원소: {self.element}')
        if not isinstance(self.position, (tuple, list)) or len(self.position) != 3:
            raise ValueError('원자마다 x, y, z 좌표 3개가 필요합니다.')
        for value in self.position:
            real(value, '원자 좌표', -10000, 10000)
        object.__setattr__(self, 'element', element)
        object.__setattr__(self, 'position', tuple(float(x) for x in self.position))


def parse_coordinates(text: str) -> tuple[Atom, ...]:
    atoms = []
    for number, line in enumerate(text.replace(';', '\n').splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 4:
            raise ValueError(f'{number}행: 원소 x y z 형식으로 입력하세요.')
        try:
            position = tuple(float(x) for x in parts[1:])
        except ValueError:
            raise ValueError(f'{number}행: 좌표는 숫자여야 합니다.') from None
        atoms.append(Atom(parts[0], position))
    if not atoms:
        raise ValueError('원자를 한 개 이상 입력하세요.')
    return tuple(atoms)


def read_xyz(path: Path) -> tuple[Atom, ...]:
    """A single standard XYZ frame, including its mandatory comment line; in Å."""
    lines = Path(path).read_text(encoding='utf-8-sig').splitlines()
    try:
        count = int(lines[0].strip())
    except (ValueError, IndexError):
        raise ValueError('XYZ 첫 행에는 원자 수가 필요합니다.') from None
    if count < 1 or len(lines) < count + 2:
        raise ValueError('XYZ 원자 수 또는 주석 행/좌표 행을 확인하세요.')
    if any(line.strip() for line in lines[count + 2:]):
        raise ValueError('한 구조의 XYZ만 지원합니다. 여러 프레임을 함께 읽을 수 없습니다.')
    atoms = parse_coordinates('\n'.join(lines[2:count + 2]))
    if len(atoms) != count:
        raise ValueError('XYZ에 선언된 원자 수와 좌표 행 수가 다릅니다.')
    return atoms


def coordinates_text(atoms):
    return '\n'.join(f'{a.element:<2} ' + ' '.join(f'{x:.12g}' for x in a.position) for a in atoms)


@dataclass(frozen=True)
class DFTConfig:
    atoms: tuple[Atom, ...] = field(default_factory=lambda: parse_coordinates(
        'O 0 0 0\nH 0 -0.757 0.587\nH 0 0.757 0.587'))
    charge: int = 0
    spin: int = 0  # N_alpha - N_beta = 2 S; NOT multiplicity
    unit: str = 'angstrom'
    basis: str = 'def2-svp'  # Default for every element not listed in basis_overrides.
    basis_overrides: dict = field(default_factory=dict)  # {element: basis}
    xc: str = 'pbe'
    dispersion: str = 'none'  # 'd4': geometry-only correction added after SCF.
    method: str = 'auto'
    scf_solver: str = 'diis'
    ecp: str = 'auto'
    relativity: str = 'auto'  # Dyall -> scalar X2C; other bases retain previous behavior.
    nuclear_model: str = 'auto'  # Dyall -> Gaussian finite nucleus; otherwise point.
    grid_scheme: str = 'auto'  # Dyall or Z>103 -> radius-independent quadrature.
    engine: str = 'pyscf'  # Existing files keep their external-engine provenance.
    calculation: str = 'single_point'  # 'gradient': also dE/dR; 'optimize': geomeTRIC minimization.
    opt_convergence: str = 'default'  # Key of OPT_CONVERGENCE.
    opt_max_steps: int = 100
    opt_constraints: str = ''  # geomeTRIC constraint text ($freeze / $set ...); ';' may separate lines.
    neb_product: str = ''  # calculation=neb: product coordinates, same atom order as atoms (Angstrom/unit).
    neb_images: int = 9  # Images including both fixed endpoints.
    neb_maxg_ev_per_angstrom: float = 0.05  # Max NEB force; the average criterion is half of it.
    scan_spec: str = ''  # calculation=scan: 'distance i j start stop points' (Angstrom / degrees, 1-based).
    freq_step_bohr: float = 0.005  # Central-difference displacement for the Hessian.
    thermo_temperature_k: float = 298.15
    thermo_pressure_pa: float = 101325.0
    symmetry_number: int = 0  # Rotational symmetry number; 0 = detect the point group.
    initial_guess: str = 'auto'
    native_density_tol: float = 1e-7
    native_residual_tol: float = 1e-6
    native_damping: float = 0.15
    native_damping_residual: float = 1e-3  # Damp only while the commutator norm exceeds this; 0 = always.
    native_stability: str = 'none'  # 'check' or 'follow' real internal orbital-Hessian stability (native only).
    native_diis_space: int = 8
    system: str = 'molecule'
    max_cycle: int = 100
    conv_tol: float = 1e-9
    grid_level: int = 3
    cube_points: int = 40
    cube_margin_bohr: float = 4.0
    threads: int = 2
    memory_mb: int = 2000

    def __post_init__(self):
        if self.system != 'molecule':
            raise ValueError('주기적 고체·표면은 PeriodicConfig와 주기계 화면을 사용하세요. 분자 엔진은 고립계만 처리합니다.')
        if self.calculation not in CALCULATIONS:
            raise ValueError('분자 계산 종류: single_point(에너지), gradient(에너지+원자 힘), '
                             'optimize(geomeTRIC 구조 최적화), frequencies(조화 진동수).')
        if self.opt_convergence not in OPT_CONVERGENCE:
            raise ValueError('구조 최적화 수렴 기준: ' + ', '.join(OPT_CONVERGENCE))
        integer(self.opt_max_steps, '구조 최적화 최대 단계', 1, 1000)
        if not isinstance(self.opt_constraints, str) or not isinstance(self.scan_spec, str):
            raise ValueError('opt_constraints와 scan_spec은 문자열이어야 합니다.')
        if self.calculation == 'scan':
            parse_scan(self.scan_spec, len(self.atoms))
        elif self.scan_spec.strip():
            raise ValueError('scan_spec은 calculation=scan에서만 사용합니다.')
        integer(self.neb_images, 'NEB 이미지 수(양 끝 포함)', 3, 41)
        real(self.neb_maxg_ev_per_angstrom, 'NEB 최대 힘 기준 [eV/Å]', 1e-4, 1)
        if not isinstance(self.neb_product, str):
            raise ValueError('neb_product는 좌표 문자열이어야 합니다.')
        if self.calculation == 'neb':
            product = parse_coordinates(self.neb_product) if self.neb_product.strip() else ()
            if [a.element for a in product] != [a.element for a in self.atoms]:
                raise ValueError('NEB 생성물(neb_product)은 반응물과 원자 수·원소 순서가 같아야 합니다.')
        elif self.neb_product.strip():
            raise ValueError('neb_product는 calculation=neb에서만 사용합니다.')
        if self.opt_constraints.strip() and self.calculation not in ('optimize', 'scan', 'transition_state'):
            raise ValueError('구조 제약은 optimize / scan / transition_state 계산에서만 사용합니다.')
        real(self.freq_step_bohr, '진동수 유한차분 변위 [bohr]', 1e-4, 0.05)
        real(self.thermo_temperature_k, '열역학 온도 [K]', 1, 5000)
        real(self.thermo_pressure_pa, '열역학 압력 [Pa]', 1, 1e9)
        integer(self.symmetry_number, '회전 대칭수 (0=자동)', 0, 1000)
        if not isinstance(self.atoms, (tuple, list)) or not 1 <= len(self.atoms) <= 100:
            raise ValueError('한 구조에 원자 1~100개를 입력하세요.')
        if not all(isinstance(atom, Atom) for atom in self.atoms):
            raise ValueError('atoms에는 원소와 3차원 좌표가 필요합니다.')
        object.__setattr__(self, 'atoms', tuple(self.atoms))
        if self.unit not in ('angstrom', 'bohr'):
            raise ValueError('좌표 단위는 angstrom 또는 bohr입니다.')
        for key, choices in (('basis', BASES), ('xc', FUNCTIONALS), ('method', ('auto', 'rks', 'uks')),
                             ('scf_solver', ('diis', 'newton')), ('ecp', ('auto', 'none', 'def2')),
                             ('relativity', ('auto', 'none', 'x2c')), ('nuclear_model', ('auto', 'point', 'gaussian')),
                             ('grid_scheme', ('auto', 'treutler', 'universal')),
                             ('engine', ('pyscf', 'native')), ('initial_guess', ('auto', 'core', 'minao', 'perturbed')),
                             ('dispersion', DISPERSIONS), ('native_stability', ('none', 'check', 'follow'))):
            value = getattr(self, key)
            if not isinstance(value, str) or value.lower().strip() not in choices:
                raise ValueError(f'지원하지 않는 {key}: {value}. 선택 가능: {", ".join(choices)}')
            object.__setattr__(self, key, value.lower().strip())
        if not isinstance(self.basis_overrides, dict):
            raise ValueError('basis_overrides는 {"원소": "기저"} 형식이어야 합니다.')
        overrides = {}
        for element, basis in self.basis_overrides.items():
            if not isinstance(element, str) or element.strip().capitalize() not in ATOMIC_NUMBERS:
                raise ValueError(f'원소별 기저: 알 수 없는 원소 {element}')
            if not isinstance(basis, str) or basis.lower().strip() not in BASES:
                raise ValueError(f'원소별 기저: 지원하지 않는 {element} 기저 {basis}. 선택 가능: {", ".join(BASES)}')
            element, basis = element.strip().capitalize(), basis.lower().strip()
            if element in overrides:
                raise ValueError(f'원소별 기저에 {element}가 두 번 지정되었습니다.')
            # An override equal to the default changes nothing physically; drop it so
            # identical calculations keep identical settings and checkpoint signatures.
            if basis != self.basis:
                overrides[element] = basis
        object.__setattr__(self, 'basis_overrides', dict(sorted(overrides.items())))
        if self.dispersion == 'd4':
            if self.xc not in D4_FUNCTIONALS:
                raise ValueError(f'D4 분산 보정의 공개 감쇠 매개변수가 {self.xc}에는 없습니다. '
                                 f'{", ".join(sorted(D4_FUNCTIONALS))} 중에서 선택하세요.')
            unsupported = sorted({a.element for a in self.atoms} - D4_ELEMENTS, key=ATOMIC_NUMBERS.get)
            if unsupported:
                raise ValueError(f'D4 참조 데이터가 없는 원소: {", ".join(unsupported)} '
                                 '(dftd4 지원: H~Lr, Cn~Og). 분산 보정을 none으로 선택하세요.')
        for element, basis in self.basis_by_element.items():
            if element in BASIS_UNAVAILABLE.get(basis, ()):
                raise ValueError(f'{basis}: {element}에 지원되는 기저 데이터가 없습니다. '
                                 f'{element} 사용 가능: '
                                 + ', '.join(available_bases((Atom(element, (0, 0, 0)),))))
        needs_ecp = any(basis in DEF2_BASES and element in DEF2_ECP_ELEMENTS
                        for element, basis in self.basis_by_element.items())
        if needs_ecp and self.ecp == 'none':
            raise ValueError('이 무거운 원소의 def2 기저는 대응하는 def2-ECP가 필요합니다. ECP를 auto 또는 def2로 선택하세요.')
        if self.ecp == 'def2' and not needs_ecp:
            raise ValueError('def2-ECP는 지원되는 무거운 원소와 def2-svp/def2-tzvp 기저를 함께 사용할 때만 선택하세요.')
        if self.resolved_ecp and self.resolved_relativity == 'x2c':
            raise ValueError(f'ECP({", ".join(self.resolved_ecp)})와 전체전자 X2C를 동시에 적용하지 않습니다. '
                             'Dyall 기저 / ECP none을 사용하세요. 원소별 기저라면 무거운 원소를 모두 Dyall로 지정하세요.')
        if self.resolved_ecp and self.resolved_nuclear_model != 'point':
            raise ValueError('ECP 계산에서는 핵 모양을 별도로 바꾸지 않습니다. nuclear_model을 auto 또는 point로 선택하세요.')
        if self.resolved_grid_scheme == 'treutler' and any(ATOMIC_NUMBERS[a.element] > 103 for a in self.atoms):
            raise ValueError('Treutler 격자의 원소별 매개변수는 Lr까지 있습니다. grid_scheme을 auto 또는 universal로 선택하세요.')
        if self.engine == 'native':
            if (any(ATOMIC_NUMBERS[a.element] > 36 for a in self.atoms) or self.resolved_ecp
                    or self.resolved_relativity != 'none' or self.resolved_nuclear_model != 'point'):
                raise ValueError('자체 엔진의 현재 범위는 H~Kr의 비상대론 전체전자·점핵입니다. '
                                 'ECP/X2C/무거운 원소는 PySCF 엔진을 선택하세요.')
            if self.scf_solver != 'diis':
                raise ValueError('자체 엔진은 직접 구현한 DIIS를 사용합니다. Newton은 PySCF 엔진에 있습니다.')
        elif self.initial_guess == 'perturbed':
            raise ValueError('perturbed 초기 추정은 자체 엔진에서만 지원합니다.')
        elif self.native_stability != 'none':
            raise ValueError('native_stability(SCF 안정성 분석)는 자체 엔진에서만 지원합니다.')
        integer(self.native_diis_space, '자체 DIIS 이력 길이', 2, 16)
        real(self.native_density_tol, '자체 밀도 수렴', 1e-12, 1e-3)
        real(self.native_residual_tol, '자체 잔차 수렴', 1e-12, 1e-3)
        real(self.native_damping, '자체 밀도 혼합 비율', 0, 0.8)
        real(self.native_damping_residual, '자체 밀도 혼합 적용 잔차', 0, 10)
        integer(self.charge, '전체 전하', -100, 100)
        integer(self.spin, '스핀 2S', 0, 200)
        if self.electrons <= 0:
            raise ValueError('전체 전하를 적용한 전자 수는 1 이상이어야 합니다.')
        if self.explicit_electrons <= 0:
            raise ValueError('ECP 코어를 제외한 명시적 전자 수는 1 이상이어야 합니다. 전하와 ECP를 확인하세요.')
        if self.spin > self.explicit_electrons or (self.explicit_electrons - self.spin) % 2:
            raise ValueError(f'계산 전자 {self.explicit_electrons}개와 스핀 2S={self.spin}이 맞지 않습니다. '
                             '계산 전자 수와 2S의 홀짝이 같아야 하고 2S는 계산 전자 수 이하여야 합니다.')
        if self.method == 'rks' and self.spin != 0:
            raise ValueError('RKS는 2S=0인 닫힌껍질만 지원합니다. auto 또는 UKS를 선택하세요.')
        for key, low, high in (('max_cycle', 1, 1000), ('grid_level', 0, 7),
                               ('cube_points', 12, 100), ('threads', 1, 64), ('memory_mb', 256, 64000)):
            integer(getattr(self, key), key, low, high)
        real(self.conv_tol, 'SCF 수렴 기준', 1e-12, 1e-3)
        real(self.cube_margin_bohr, '밀도 상자 여백 [bohr]', 1, 20)
        scale = 0.529177210903 if self.unit == 'bohr' else 1.0
        for i, atom in enumerate(self.atoms):
            for other in self.atoms[:i]:
                if math.dist(atom.position, other.position) * scale < 0.1:
                    raise ValueError('원자 간 거리가 0.1 Å 미만입니다. 중복 좌표/단위를 확인하세요.')

    @property
    def electrons(self):
        return sum(ATOMIC_NUMBERS[a.element] for a in self.atoms) - self.charge

    @property
    def basis_by_element(self):
        """The basis actually used for each element present in the structure."""
        return {element: self.basis_overrides.get(element, self.basis)
                for element in sorted({a.element for a in self.atoms})}

    @property
    def basis_label(self):
        return describe_bases(self.basis, self.basis_by_element)

    @property
    def unused_basis_overrides(self):
        return sorted(set(self.basis_overrides) - {a.element for a in self.atoms})

    @property
    def uses_relativistic_basis(self):
        return any(basis in RELATIVISTIC_BASES for basis in self.basis_by_element.values())

    @property
    def resolved_ecp(self):
        if self.ecp == 'none':
            return {}
        return {symbol: basis for symbol, basis in self.basis_by_element.items()
                if basis in DEF2_BASES and symbol in DEF2_ECP_ELEMENTS}

    @property
    def ecp_core_electrons(self):
        mapping = self.resolved_ecp
        return sum(DEF2_CORE_ELECTRONS[atom.element] for atom in self.atoms if atom.element in mapping)

    @property
    def resolved_relativity(self):
        return ('x2c' if self.uses_relativistic_basis else 'none') if self.relativity == 'auto' else self.relativity

    @property
    def resolved_nuclear_model(self):
        return ('gaussian' if self.uses_relativistic_basis else 'point') if self.nuclear_model == 'auto' else self.nuclear_model

    @property
    def resolved_grid_scheme(self):
        if self.grid_scheme != 'auto':
            return self.grid_scheme
        return ('universal' if self.uses_relativistic_basis
                or any(ATOMIC_NUMBERS[a.element] > 103 for a in self.atoms) else 'treutler')

    @property
    def explicit_electrons(self):
        return self.electrons - self.ecp_core_electrons

    @property
    def resolved_method(self):
        return ('rks' if self.spin == 0 else 'uks') if self.method == 'auto' else self.method

    @property
    def formula(self):
        counts = {}
        for atom in self.atoms:
            counts[atom.element] = counts.get(atom.element, 0) + 1
        order = (['C'] + (['H'] if 'H' in counts else [])
                 + sorted(set(counts) - {'C', 'H'})) if 'C' in counts else sorted(counts)
        return ''.join(e + (str(counts[e]) if counts[e] > 1 else '') for e in order)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict):
            raise ValueError('분자 설정 파일에는 JSON 객체가 필요합니다.')
        unknown = set(data) - set(cls.__dataclass_fields__)
        if unknown:
            raise ValueError('알 수 없는 분자 설정: ' + ', '.join(sorted(unknown)))
        if 'atoms' not in data:
            raise ValueError('분자 설정에 atoms가 필요합니다. 수소 검증용 설정과 구분하세요.')
        data = dict(data)
        try:
            data['atoms'] = tuple(Atom(**item) for item in data['atoms'])
            return cls(**data)
        except (TypeError, KeyError) as exc:
            raise ValueError(f'분자 설정 형식 오류: {exc}') from None


def load_dft_config(path):
    return DFTConfig.from_dict(json.loads(Path(path).read_text(encoding='utf-8-sig')))


def save_dft_config(config, path):
    Path(path).write_text(json.dumps(config.to_dict(), ensure_ascii=False, indent=2, allow_nan=False),
                          encoding='utf-8')
