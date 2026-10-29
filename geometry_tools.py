"""Engine-independent geometry measurements and explicit numerical presets."""
from dataclasses import replace
import numpy as np

BOHR_ANGSTROM = 0.529177210903
PRESETS = {
    '빠른 확인': dict(grid_level=1, conv_tol=1e-6, max_cycle=60, cube_points=24,
                    native_density_tol=1e-5, native_residual_tol=1e-4),
    '정밀 격자': dict(grid_level=4, conv_tol=1e-10, max_cycle=200, cube_points=60,
                    native_density_tol=1e-8, native_residual_tol=1e-7),
}


def apply_preset(config, name):
    return replace(config, **PRESETS[name])


def measure(atoms, indices, unit='angstrom'):
    if unit not in ('angstrom', 'bohr') or len(indices) not in (2, 3):
        raise ValueError('원자 번호2개는 거리, 3개는 가운데 원자의 각도입니다. 단위는 angstrom 또는 bohr입니다.')
    if any(not isinstance(i, int) or not 1 <= i <= len(atoms) for i in indices) or len(set(indices)) != len(indices):
        raise ValueError('서로 다른 원자 번호를1부터 원자 수까지 입력하세요.')
    p = np.array([atoms[i - 1].position for i in indices]) * (BOHR_ANGSTROM if unit == 'bohr' else 1.)
    if len(p) == 2:
        return float(np.linalg.norm(p[0] - p[1])), 'angstrom'
    a, b = p[0] - p[1], p[2] - p[1]
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    if norm < 1e-14:
        raise ValueError('겹친 원자의 각도는 정의되지 않습니다.')
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b) / norm, -1, 1)))), 'degree'
