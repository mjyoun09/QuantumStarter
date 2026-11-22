"""Periodic PBE with GTH pseudopotentials and 2D/3D Coulomb boundaries."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import physical_constants
from molecular import ATOMIC_NUMBERS, save_dft_config
from workflow import VERSION

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


BOHR_ANGSTROM = physical_constants['Bohr radius'][0] * 1e10
HARTREE_PER_BOHR3_GPA = (physical_constants['Hartree energy'][0]
                         / physical_constants['Bohr radius'][0]**3 / 1e9)


def periodic_gradient(mf, config, cell):
    """Forces and stress per input cell from PySCF's FFTDF derivative drivers.

    stress = (1/V) dX/d(strain) with X = E, or the Mermin free energy F = E - TS with
    smearing; pressure = -trace/3. Both were checked against finite differences.
    """
    solver = mf.nuc_grad_method()
    gradient = np.asarray(solver.kernel())
    stress = np.asarray(solver.get_stress())
    if not (np.isfinite(gradient).all() and np.isfinite(stress).all()):
        raise RuntimeError('유한한 주기계 기울기/응력을 얻지 못했습니다.')
    forces = -gradient * HARTREE_TO_EV / BOHR_ANGSTROM
    stress_gpa = stress * HARTREE_PER_BOHR3_GPA
    return {'driver': 'PySCF pbc analytic nuclear gradient / get_stress (FFTDF, uniform XC grid)',
            'stress_method_note': 'PySCF 2.14 stress includes finite differences of kinetic and overlap integrals (strain 1e-5).',
            'derivative_of': 'Mermin free energy F = E - TS' if config.smearing_ev else 'total energy E',
            'gradient_hartree_per_bohr': gradient.tolist(), 'forces_ev_per_angstrom': forces.tolist(),
            'max_atomic_force_ev_per_angstrom': float(np.linalg.norm(forces, axis=1).max()),
            'net_force_ev_per_angstrom': forces.sum(axis=0).tolist(),
            'stress_hartree_per_bohr3': stress.tolist(), 'stress_gpa': stress_gpa.tolist(),
            'pressure_gpa': float(-np.trace(stress_gpa) / 3),
            'stress_convention': 'sigma = (1/V) dX/d(strain); positive diagonal = energy rises on expansion',
            'fft_mesh': [int(n) for n in cell.mesh], 'cell_volume_bohr3': float(cell.vol)}


def calculate_periodic(config, directory, progress=print, cancel_directory=None):
    import pyscf
    from pyscf import lib
    from pyscf.pbc import dft, scf
    from scf_control import check_cancelled
    started = time.perf_counter()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    cancel_directory = directory if cancel_directory is None else Path(cancel_directory)
    check_cancelled(cancel_directory)
    save_dft_config(config, directory / 'settings.json')
    lib.num_threads(config.threads)
    cell, valence = build_cell(config)
    kpts = cell.make_kpts(config.kmesh)
    nk = len(kpts)
    unrestricted = config.spin != 0
    fftdf = config.integral_scheme == 'fftdf'
    if fftdf and config.fft_mesh:
        cell.mesh = np.asarray(config.fft_mesh)
        cell.build(False, False)
    mf = dft.KUKS(cell, kpts) if unrestricted else dft.KRKS(cell, kpts)
    if not fftdf:
        mf = mf.density_fit()
    if unrestricted:
        # PySCF's default spin convention is over the k-point supercell.
        mf.nelec = ((cell.nelectron + config.spin) * nk // 2,
                    (cell.nelectron - config.spin) * nk // 2)
    mf.xc = config.xc
    if fftdf:
        # Coulomb and XC on the same uniform FFT mesh; energy finite differences
        # check both gradient and stress (VALIDATION 0.9.1/0.9.2).
        mf.grids = dft.gen_grid.UniformGrids(cell)
        mf.grids.mesh = np.asarray(cell.mesh)
        mf.small_rho_cutoff = 0
    else:
        mf.grids = dft.gen_grid.BeckeGrids(cell)
        mf.grids.level = config.grid_level
    if config.smearing_ev:
        mf = scf.addons.smearing_(mf, sigma=config.smearing_ev / HARTREE_TO_EV, method='fermi', fix_spin=unrestricted)
    mf.max_cycle, mf.conv_tol = config.max_cycle, config.conv_tol
    mf.conv_tol_grad = np.sqrt(config.conv_tol)
    mf.chkfile = str(directory / 'scf.chk')
    history = []

    def callback(env):
        check_cancelled(cancel_directory)
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
    gradient = None
    check_cancelled(cancel_directory)
    if config.calculation == 'gradient':
        if mf.converged:
            progress('주기계 힘·응력(PySCF FFTDF 해석적 기울기)을 계산합니다.')
            gradient = periodic_gradient(mf, config, cell)
        else:
            progress('SCF 미수렴: 주기계 힘·응력을 계산하지 않았습니다.')
    check_cancelled(cancel_directory)
    from .periodic_analysis import band_structure, density_of_states
    dos = density_of_states(mf.mo_energy, mf.mo_occ, [1 / nk] * nk, unrestricted, config.dos_sigma_ev, HARTREE_TO_EV)
    bands = None
    if config.band_path and mf.converged:
        progress('밴드 구조: 수렴한 밀도로 고대칭 경로의 고유값을 계산합니다(비자기일관).')
        bands = band_structure(mf, cell, config, unrestricted, HARTREE_TO_EV)
    check_cancelled(cancel_directory)
    warnings = ['SCF 수렴은 기저·k점·수치 격자·진공·슬래브 두께 수렴의 보증이 아닙니다.',
                'GTH 의사퍼텐셜 계산이며 SOC·분산 보정은 포함하지 않습니다.']
    if not fftdf:
        warnings.append('GDF 경로에서는 주기계 힘·응력·셀 최적화를 지원하지 않습니다.')
    if config.smearing_ev:
        warnings.append('내부에너지 E, 자유에너지 F=E−TS, 외삽 에너지 (E+F)/2를 구분하세요. 비교에는 같은 smearing을 사용하세요.')
    if not mf.converged:
        warnings.append('SCF 미수렴: 후보 선정이나 정량 비교에 사용하지 마세요.')
    if config.boundary == 'slab':
        warnings.append('xy 주기·z 비주기 2D Coulomb 경계입니다. z 방향 셀 길이에 대한 수치 수렴도 확인하세요.')
    result = {
        'schema_version': 3, 'version': VERSION, 'engine': 'pyscf-pbc', 'pyscf_version': pyscf.__version__,
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
        'integral_scheme': config.integral_scheme.upper(), 'fft_mesh': [int(n) for n in cell.mesh] if fftdf else None,
        'xc_grid': 'uniform FFT mesh' if fftdf else f'Becke level {config.grid_level}',
        'calculation': config.calculation, 'gradient': gradient, 'coulomb_dimension': cell.dimension, 'scf_history': history, 'warnings': warnings,
        'restarted_same_settings': restarted, 'dos': dos, 'band_structure': bands,
    }
    with (directory / 'dos.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['energy_ev'] + [f'dos_spin{s}' for s in range(len(dos['dos_by_spin']))])
        writer.writerows([e, *values] for e, *values in zip(dos['energy_grid_ev'], *dos['dos_by_spin']))
    if bands:
        with (directory / 'bands.csv').open('w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['k_index', 'axis_inverse_angstrom', 'spin', 'band', 'energy_ev'])
            for s, channel in enumerate(bands['bands_ev_by_spin']):
                for k, (x, row) in enumerate(zip(bands['axis_inverse_angstrom'], channel)):
                    writer.writerows([k, x, s, n, e] for n, e in enumerate(row))
    np.savez_compressed(directory / 'electronic_structure.npz', density_matrix=dm, overlap=overlap,
                        mo_coeff=np.asarray(mf.mo_coeff), mo_energy_hartree=np.asarray(mf.mo_energy),
                        mo_occupation=np.asarray(mf.mo_occ), kpoints_inverse_bohr=kpts, lattice_angstrom=config.lattice)
    with (directory / 'scf.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['cycle', 'energy_hartree', 'gradient_norm'])
        writer.writeheader()
        writer.writerows(history)
    (directory / 'results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    return result
