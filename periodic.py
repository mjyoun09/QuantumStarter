"""Periodic inputs: lattice rows and Cartesian positions in angstrom, spin per cell."""
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import numpy as np
from molecular import ATOMIC_NUMBERS, integer, real
from element_catalog import periodic_supported_elements

PERIODIC_BASES = ('gth-dzvp-molopt-sr', 'gth-szv-molopt-sr', 'gth-dzvp', 'gth-tzvp')


@dataclass(frozen=True)
class PeriodicAtom:
    element: str
    position: tuple[float, float, float]

    def __post_init__(self):
        if not isinstance(self.element, str) or self.element.strip().capitalize() not in ATOMIC_NUMBERS:
            raise ValueError(f'알 수 없는 원소: {self.element}')
        if not isinstance(self.position, (list, tuple)) or len(self.position) != 3:
            raise ValueError('각 원자의 Cartesian 좌표 3개가 필요합니다.')
        for x in self.position:
            real(x, '좌표 (Å)', -10000, 10000)
        object.__setattr__(self, 'element', self.element.strip().capitalize())
        object.__setattr__(self, 'position', tuple(float(x) for x in self.position))


@dataclass(frozen=True)
class PeriodicConfig:
    atoms: tuple[PeriodicAtom, ...]
    lattice: tuple[tuple[float, float, float], ...]
    system: str = 'periodic'
    boundary: str = 'bulk'
    kmesh: tuple[int, int, int] = (1, 1, 1)
    charge: int = 0
    spin: int = 0  # fixed Nalpha-Nbeta PER INPUT CELL, including smearing
    basis: str = 'gth-dzvp-molopt-sr'
    pseudo: str = 'gth-pbe'
    xc: str = 'pbe'
    smearing_ev: float = 0.1
    precision: float = 1e-7
    conv_tol: float = 1e-8
    max_cycle: int = 100
    grid_level: int = 3
    threads: int = 4
    memory_mb: int = 4000
    dos_sigma_ev: float = 0.1  # Gaussian broadening of the k-mesh DOS.
    band_path: str = ''  # '' = no band structure; 'auto' = ASE default path; or e.g. 'GXWKGL'.
    band_points: int = 60
    integral_scheme: str = 'gdf'  # 'gdf' (Becke XC grid; slabs) or 'fftdf' (bulk; uniform FFT XC grid).
    fft_mesh: tuple = ()  # FFTDF only: explicit FFT mesh; () = PySCF's choice from `precision`.
    calculation: str = 'single_point'  # gradient / cell_optimize: FFTDF bulk only.
    cell_relaxation: str = 'full'  # full cell shape, or isotropic volume with fixed shape.
    cell_max_steps: int = 50
    cell_force_tol: float = 0.02  # eV/angstrom, maximum atomic force norm.
    cell_stress_tol_gpa: float = 0.1
    cell_energy_tol_ev: float = 1e-5  # E or Mermin F change per input cell.

    def __post_init__(self):
        if self.system != 'periodic':
            raise ValueError('PeriodicConfig의 system은 periodic이어야 합니다.')
        if not isinstance(self.atoms, (tuple, list)) or not 1 <= len(self.atoms) <= 200:
            raise ValueError('주기계 입력은 원자 1~200개가 필요합니다.')
        if not all(isinstance(a, PeriodicAtom) for a in self.atoms):
            raise ValueError('주기계 원자 형식을 확인하세요.')
        object.__setattr__(self, 'atoms', tuple(self.atoms))
        cell = np.asarray(self.lattice, dtype=float)
        if cell.shape != (3, 3) or not np.isfinite(cell).all() or np.max(np.abs(cell)) > 1000:
            raise ValueError('lattice는 유한한 3×3 격자 벡터(각 행, Å)여야 합니다.')
        if np.linalg.det(cell) <= 1e-6 or np.linalg.cond(cell) > 1e5:
            raise ValueError('격자는 오른손 좌표계이며 부피가 양수여야 합니다.')
        object.__setattr__(self, 'lattice', tuple(tuple(float(x) for x in row) for row in cell))
        if self.boundary not in ('bulk', 'slab'):
            raise ValueError('boundary는 bulk 또는 slab입니다.')
        if not isinstance(self.kmesh, (list, tuple)) or len(self.kmesh) != 3:
            raise ValueError('kmesh에 정수 3개가 필요합니다.')
        for k in self.kmesh:
            integer(k, 'k점 격자', 1, 12)
        if np.prod(self.kmesh) > 216:
            raise ValueError('현재 최대 k점 수는 216입니다.')
        object.__setattr__(self, 'kmesh', tuple(self.kmesh))
        if self.boundary == 'slab':
            if not np.allclose(cell[:2, 2], 0, atol=1e-10) or not np.allclose(cell[2, :2], 0, atol=1e-10):
                raise ValueError('slab은 xy 평면, 세 번째 격자 벡터는 z 방향이어야 합니다.')
            if self.kmesh[2] != 1:
                raise ValueError('slab의 z 방향 k점은 1이어야 합니다.')
            z = np.asarray([a.position[2] for a in self.atoms])
            if z.min() < 0 or z.max() > cell[2, 2] or cell[2, 2] - np.ptp(z) < 6:
                raise ValueError('slab 원자는 셀 내부에 두고 총 진공 간격을 6 Å 이상 확보하세요.')
        integer(self.charge, '전체 전하', 0, 0)
        integer(self.spin, '셀당 스핀 2S', 0, 200)
        if self.basis not in PERIODIC_BASES or self.pseudo != 'gth-pbe' or self.xc != 'pbe':
            raise ValueError('현재 주기계는 PBE / gth-pbe / 지정된 GTH 기저 조합을 지원합니다.')
        missing = sorted({a.element for a in self.atoms} - periodic_supported_elements(self.basis))
        if missing:
            raise ValueError(f'{", ".join(missing)}: {self.basis} / {self.pseudo} 주기계 데이터가 없습니다. '
                             '원소 입력 지원과 주기계 계산 지원은 다릅니다. 고립계는 분자 화면의 Dyall 기저를 사용할 수 있습니다.')
        for key, low, high in (('smearing_ev', 0, 1), ('precision', 1e-12, 1e-4), ('conv_tol', 1e-12, 1e-3)):
            real(getattr(self, key), key, low, high)
        for key, low, high in (('max_cycle', 1, 1000), ('grid_level', 0, 6), ('threads', 1, 64), ('memory_mb', 256, 64000)):
            integer(getattr(self, key), key, low, high)
        real(self.dos_sigma_ev, 'DOS 확장 폭 [eV]', 1e-3, 2)
        if self.integral_scheme not in ('gdf', 'fftdf'):
            raise ValueError('integral_scheme은 gdf 또는 fftdf입니다.')
        if self.calculation not in ('single_point', 'gradient', 'cell_optimize'):
            raise ValueError('주기계 calculation은 single_point, gradient, cell_optimize입니다.')
        if self.integral_scheme == 'fftdf' and self.boundary != 'bulk':
            raise ValueError('FFTDF 경로는 bulk 전용입니다. slab의 2D Coulomb 처리는 GDF 경로에만 있습니다.')
        if self.calculation in ('gradient', 'cell_optimize') and self.integral_scheme != 'fftdf':
            raise ValueError('주기계 힘·응력은 integral_scheme=fftdf(bulk)에서만 계산합니다. '
                             'PySCF 2.14의 GDF는 주기계 기울기를 지원하지 않습니다.')
        if not isinstance(self.fft_mesh, (tuple, list)) or len(self.fft_mesh) not in (0, 3):
            raise ValueError('fft_mesh는 비우거나 정수 3개입니다.')
        for n in self.fft_mesh:
            integer(n, 'FFT mesh', 5, 400)
        if self.fft_mesh and self.integral_scheme != 'fftdf':
            raise ValueError('fft_mesh는 integral_scheme=fftdf에서만 씁니다.')
        object.__setattr__(self, 'fft_mesh', tuple(self.fft_mesh))
        if self.cell_relaxation not in ('full', 'isotropic'):
            raise ValueError('cell_relaxation은 full 또는 isotropic입니다.')
        integer(self.cell_max_steps, '셀 최적화 최대 단계', 1, 500)
        real(self.cell_force_tol, '원자 힘 기준 [eV/Å]', 1e-6, 1)
        real(self.cell_stress_tol_gpa, '응력 기준 [GPa]', 1e-5, 10)
        real(self.cell_energy_tol_ev, '셀 에너지 변화 기준 [eV]', 1e-10, .1)
        integer(self.band_points, '밴드 경로 점 수', 10, 500)
        if not isinstance(self.band_path, str) or not all(c.isalnum() or c == ',' for c in self.band_path):
            raise ValueError("band_path는 ''(없음), 'auto', 또는 GXWKGL 같은 고대칭점 문자열입니다.")
        from ase.geometry import find_mic
        pos = np.asarray([a.position for a in self.atoms])
        if len(pos) > 1:
            delta = pos[:, None, :] - pos[None, :, :]
            _, distances = find_mic(delta.reshape(-1, 3), cell, pbc=self.pbc)
            distances = distances.reshape(len(pos), len(pos))
            np.fill_diagonal(distances, np.inf)
            if distances.min() < .1:
                raise ValueError('주기 경계에서 중복되거나 지나치게 가까운 원자가 있습니다.')

    @property
    def pbc(self):
        return (True, True, self.boundary == 'bulk')

    @property
    def formula(self):
        from collections import Counter
        counts = Counter(a.element for a in self.atoms)
        return ''.join(e + (str(n) if n > 1 else '') for e, n in sorted(counts.items()))

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict) or 'atoms' not in data or 'lattice' not in data:
            raise ValueError('주기계 설정에는 atoms와 lattice가 필요합니다.')
        if set(data) - set(cls.__dataclass_fields__):
            raise ValueError('알 수 없는 주기계 설정: ' + ', '.join(sorted(set(data) - set(cls.__dataclass_fields__))))
        try:
            data = dict(data)
            data['atoms'] = tuple(PeriodicAtom(**a) for a in data['atoms'])
            return cls(**data)
        except (TypeError, KeyError) as exc:
            raise ValueError(f'주기계 설정 형식 오류: {exc}') from exc

    @classmethod
    def from_ase(cls, atoms, **settings):
        pbc = tuple(bool(x) for x in atoms.pbc)
        boundary = settings.pop('boundary', 'slab' if pbc == (True, True, False) else 'bulk')
        if pbc not in ((True, True, True), (True, True, False)):
            raise ValueError('구조의 주기 경계는 xyz 또는 xy여야 합니다.')
        return cls(atoms=tuple(PeriodicAtom(e, tuple(p)) for e, p in zip(atoms.symbols, atoms.positions)),
                   lattice=tuple(map(tuple, atoms.cell.array)), boundary=boundary, **settings)


def config_from_dict(data):
    from molecular import DFTConfig
    if not isinstance(data, dict):
        raise ValueError('설정에는 JSON 객체가 필요합니다.')
    return (PeriodicConfig if data.get('system') == 'periodic' else DFTConfig).from_dict(data)


def load_calculation_config(path):
    return config_from_dict(json.loads(Path(path).read_text(encoding='utf-8-sig')))


def read_periodic_structure(path, **settings):
    from ase.io import read
    path = Path(path)
    fmt = 'cif' if path.suffix.lower() == '.cif' else 'extxyz' if path.suffix.lower() in ('.xyz', '.extxyz') else 'vasp'
    structures = read(str(path), index=':', format=fmt)
    if len(structures) != 1:
        raise ValueError('한 구조만 포함한 POSCAR / CIF / extxyz가 필요합니다.')
    return PeriodicConfig.from_ase(structures[0], **settings)
