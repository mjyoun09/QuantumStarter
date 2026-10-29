"""Portable, non-pickle SCF checkpoints and cooperative cancellation."""
import hashlib
import json
from pathlib import Path
import os
import numpy as np


class CalculationCancelled(RuntimeError):
    pass


def check_cancelled(directory):
    if (Path(directory) / 'cancel.request').exists():
        raise CalculationCancelled('사용자가 계산 중단을 요청했습니다. 마지막 저장 체크포인트를 확인하세요.')


def signature(config):
    # Numerical physics and basis order must match. Resources, stopping criteria
    # and SCF iteration strategy may be changed on restart. Dispersion is a
    # geometry-only correction added after SCF; it never changes the density.
    data = config.to_dict()
    for key in ('engine', 'calculation', 'max_cycle', 'conv_tol', 'threads', 'memory_mb', 'cube_points', 'cube_margin_bohr',
                'initial_guess', 'scf_solver', 'native_density_tol', 'native_residual_tol',
                'native_damping', 'native_diis_space', 'dispersion', 'opt_convergence', 'opt_max_steps',
                'freq_step_bohr', 'thermo_temperature_k'):
        data.pop(key, None)
    # Only overrides for elements in the structure change the AO basis. Without any,
    # the hash equals the one written before per-element bases existed.
    used = {e: b for e, b in data.pop('basis_overrides', {}).items() if e in config.basis_by_element}
    if used:
        data['basis_overrides'] = used
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def save_checkpoint(directory, config, dm, cycle):
    path = Path(directory) / 'checkpoint.npz'
    temporary = path.with_suffix('.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, density_matrix=np.asarray(dm), signature=np.array(signature(config)),
                            cycle=np.array(int(cycle)), schema_version=np.array(2), source_engine=np.array(config.engine))
    os.replace(temporary, path)


def load_checkpoint(path, config, overlap):
    try:
        with np.load(path, allow_pickle=False) as data:
            if int(data['schema_version']) != 2 or str(data['signature']) != signature(config):
                raise ValueError('체크포인트 버전 또는 구조·기저·XC·전하/스핀·격자가 현재 설정과 다릅니다.')
            dm = np.asarray(data['density_matrix'], dtype=float)
        n = overlap.shape[0]
        expected = (n, n) if config.resolved_method == 'rks' else (2, n, n)
        if dm.shape != expected or not np.isfinite(dm).all() or not np.allclose(dm, dm.swapaxes(-1, -2), atol=1e-9):
            raise ValueError('체크포인트 밀도행렬의 크기·유한성·대칭성을 확인하세요.')
        matrices = [dm] if dm.ndim == 2 else dm
        targets = [config.explicit_electrons] if dm.ndim == 2 else [
            (config.explicit_electrons + config.spin) / 2, (config.explicit_electrons - config.spin) / 2]
        se, su = np.linalg.eigh(overlap)
        root_s = (su * np.sqrt(se)) @ su.T
        for matrix, target in zip(matrices, targets):
            if abs(np.einsum('ij,ji->', matrix, overlap) - target) > 1e-6:
                raise ValueError('체크포인트의 스핀별 전자 수가 설정과 맞지 않습니다.')
            occupation = np.linalg.eigvalsh(root_s @ matrix @ root_s)
            if occupation.min() < -1e-7 or occupation.max() > (2 if dm.ndim == 2 else 1) + 1e-7:
                raise ValueError('체크포인트 밀도행렬의 점유수 범위가 잘못되었습니다.')
        return dm
    except (OSError, KeyError) as exc:
        raise ValueError(f'체크포인트를 읽을 수 없습니다: {exc}') from exc
