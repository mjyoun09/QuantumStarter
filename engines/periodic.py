"""Periodic PBE with GTH pseudopotentials and 2D/3D Coulomb boundaries."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import physical_constants
from molecular import ATOMIC_NUMBERS, save_dft_config

HARTREE_TO_EV = physical_constants['Hartree energy in eV'][0]


def build_cell(config):
    from pyscf.pbc import gto
    from pyscf.pbc.gto import basis, pseudo
    elements = sorted({a.element for a in config.atoms})
    bases, potentials = {}, {}
    for element in elements:
        try:
            bases[element] = basis.load(config.basis, element)
            potentials[element] = pseudo.load(config.pseudo, element)
            if not bases[element] or not potentials[element]:
                raise ValueError('empty data')
        except Exception as exc:
            raise ValueError(f'{element}: {config.basis} / {config.pseudo} 데이터가 없습니다.') from exc
    valence = {e: sum(p[0]) for e, p in potentials.items()}
    ne = sum(valence[a.element] for a in config.atoms) - config.charge
    if config.spin > ne or (ne - config.spin) % 2:
        raise ValueError(f'명시적 전자 {ne}개와 셀당 2S={config.spin}의 홀짝/크기가 맞지 않습니다.')
    cell = gto.Cell()
    cell.atom = [(a.element, a.position) for a in config.atoms]
    cell.a, cell.unit = config.lattice, 'angstrom'
    cell.basis, cell.pseudo = bases, potentials
    cell.charge, cell.spin = config.charge, config.spin
    cell.dimension = 2 if config.boundary == 'slab' else 3
    cell.low_dim_ft_type = None
    cell.precision, cell.max_memory, cell.verbose = config.precision, config.memory_mb, 0
    cell.build()
    if cell.nelectron != ne:
        raise RuntimeError('의사퍼텐셜의 원자가 전자 수가 엔진과 일치하지 않습니다.')
    if cell.nao_nr() > 2000:
        raise ValueError('현재 주기계 실행 한도는 기저함수 2000개입니다.')
    if (ne + config.spin) // 2 > cell.nao_nr():
        raise ValueError('기저함수에 전자를 배치할 공간이 부족합니다.')
    return cell, valence


def calculate_periodic(config, directory, progress=print):
    import pyscf
    from pyscf import lib
    from pyscf.pbc import dft, scf
    started = time.perf_counter()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    save_dft_config(config, directory / 'settings.json')
    lib.num_threads(config.threads)
    cell, valence = build_cell(config)
    kpts = cell.make_kpts(config.kmesh)
    nk = len(kpts)
    unrestricted = config.spin != 0
    mf = (dft.KUKS(cell, kpts) if unrestricted else dft.KRKS(cell, kpts)).density_fit()
    if unrestricted:
        # PySCF's default spin convention is over the k-point supercell.
        mf.nelec = ((cell.nelectron + config.spin) * nk // 2,
                    (cell.nelectron - config.spin) * nk // 2)
    mf.xc = config.xc
    mf.grids = dft.gen_grid.BeckeGrids(cell)
    mf.grids.level = config.grid_level
    if config.smearing_ev:
        mf = scf.addons.smearing_(mf, sigma=config.smearing_ev / HARTREE_TO_EV, method='fermi', fix_spin=unrestricted)
    mf.max_cycle, mf.conv_tol = config.max_cycle, config.conv_tol
    mf.conv_tol_grad = np.sqrt(config.conv_tol)
    mf.chkfile = str(directory / 'scf.chk')
    history = []

    def callback(env):
        row = {'cycle': len(history) + 1, 'energy_hartree': float(env['e_tot']),
               'gradient_norm': float(env.get('norm_gorb', 0))}
        history.append(row)
        progress(f'SCF {row["cycle"]:3d} · E = {row["energy_hartree"]:.12f} Ha / cell')

    mf.callback = callback
    progress(f'{config.formula} · {config.boundary} · k={config.kmesh} · {cell.nao_nr()} AO · '
             f'명시적 전자 {cell.nelectron}개 / cell · PBE / {config.basis}')
    energy = float(mf.kernel())
    # PySCF may fail its final extra-cycle check after early stopping. One restart
    # from the final density uses identical tolerances and only the remaining budget.
    restarted = False
    if not mf.converged and 0 < len(history) < config.max_cycle:
        restarted = True
        mf.max_cycle = config.max_cycle - len(history)
        progress('최종 수렴 확인 실패: 남은 반복 한도에서 동일 조건·최종 밀도로 한 번 재시도합니다.')
        energy = float(mf.kernel(dm0=mf.make_rdm1()))
    if not np.isfinite(energy):
        raise RuntimeError('유한한 주기계 에너지를 얻지 못했습니다.')
    free_energy = float(mf.e_free) if config.smearing_ev else energy
    extrapolated = (energy + free_energy) / 2
    dm = np.asarray(mf.make_rdm1())
    overlap = np.asarray(mf.get_ovlp())
    total_dm = dm.sum(axis=0) if unrestricted else dm
    ne_trace = float(np.einsum('kij,kji->', total_dm, overlap).real / nk)
    spin_trace = float(np.einsum('kij,kji->', dm[0] - dm[1], overlap).real / nk) if unrestricted else 0.
    if abs(ne_trace - cell.nelectron) > 1e-5 or abs(spin_trace - config.spin) > 1e-5:
        raise RuntimeError(f'전자 수/스핀 검증 실패: Tr(DS)={ne_trace}, spin={spin_trace}')
    warnings = ['SCF 수렴은 기저·k점·수치 격자·진공·슬래브 두께 수렴의 보증이 아닙니다.',
                '고정 구조의 단일점 에너지입니다. 현재 주기계 DFT는 힘·응력·구조 최적화를 출력하지 않습니다.',
                'GTH 의사퍼텐셜 계산이며 SOC·분산 보정은 포함하지 않습니다.']
    if config.smearing_ev:
        warnings.append('내부에너지 E, 자유에너지 F=E−TS, 외삽 에너지 (E+F)/2를 구분하세요. 비교에는 같은 smearing을 사용하세요.')
    if not mf.converged:
        warnings.append('SCF 미수렴: 후보 선정이나 정량 비교에 사용하지 마세요.')
    if config.boundary == 'slab':
        warnings.append('xy 주기·z 비주기 2D Coulomb 경계입니다. z 방향 셀 길이에 대한 수치 수렴도 확인하세요.')
    result = {
        'schema_version': 3, 'version': '0.7.0', 'engine': 'pyscf-pbc', 'pyscf_version': pyscf.__version__,
        'created_at_utc': datetime.now(timezone.utc).isoformat(), 'config': config.to_dict(),
        'formula': config.formula, 'status': 'converged' if mf.converged else 'not_converged',
        'converged': bool(mf.converged), 'total_energy_hartree': energy, 'total_energy_ev': energy * HARTREE_TO_EV,
        'free_energy_hartree': free_energy, 'free_energy_ev': free_energy * HARTREE_TO_EV,
        'extrapolated_energy_hartree': extrapolated, 'extrapolated_energy_ev': extrapolated * HARTREE_TO_EV,
        'energy_normalization': 'per_input_cell', 'explicit_electrons': cell.nelectron,
        'electrons': sum(ATOMIC_NUMBERS[a.element] for a in config.atoms) - config.charge,
        'valence_electrons_by_element': valence, 'electron_trace': ne_trace, 'spin_trace_per_cell': spin_trace,
        'kpoints_scaled': cell.get_scaled_kpts(kpts).tolist(), 'kpoint_weights': [1 / nk] * nk,
        'nao': cell.nao_nr(), 'iterations': len(history), 'elapsed_seconds': time.perf_counter() - started,
        'integral_scheme': 'GDF', 'coulomb_dimension': cell.dimension, 'scf_history': history, 'warnings': warnings,
        'restarted_same_settings': restarted,
    }
    np.savez_compressed(directory / 'electronic_structure.npz', density_matrix=dm, overlap=overlap,
                        mo_coeff=np.asarray(mf.mo_coeff), mo_energy_hartree=np.asarray(mf.mo_energy),
                        mo_occupation=np.asarray(mf.mo_occ), kpoints_inverse_bohr=kpts, lattice_angstrom=config.lattice)
    with (directory / 'scf.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['cycle', 'energy_hartree', 'gradient_norm'])
        writer.writeheader()
        writer.writerows(history)
    (directory / 'results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    return result
