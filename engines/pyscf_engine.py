"""Common molecular setup/export with explicit native or external PySCF SCF."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import time

import numpy as np
from scipy.constants import physical_constants

from molecular import ATOMIC_NUMBERS, DEF2_CORE_ELECTRONS, TRANSITION_METALS, save_dft_config
from scf_control import CalculationCancelled, check_cancelled, save_checkpoint, load_checkpoint

HARTREE_TO_EV = physical_constants['Hartree energy in eV'][0]


def configure_grids(grids, config):
    """Avoid incomplete element-radius tables for the all-element basis path."""
    from pyscf.dft import gen_grid, radi
    grids.level = config.grid_level
    grids.becke_scheme = gen_grid.original_becke
    if config.resolved_grid_scheme == 'universal':
        grids.radi_method = radi.gauss_chebyshev
        grids.prune = None
        grids.radii_adjust = None
        grids.atomic_radii = None
        # Heavy all-electron densities need more radial resolution than the
        # standard period-based grid. Keep the level sequence increasing and
        # record the actual counts; this is not an energy-convergence guarantee.
        grids.atom_grid = {a.element: (3 * gen_grid._default_rad(ATOMIC_NUMBERS[a.element], config.grid_level),
                                       gen_grid._default_ang(ATOMIC_NUMBERS[a.element], config.grid_level))
                           for a in config.atoms if ATOMIC_NUMBERS[a.element] > 36}
    else:
        grids.radi_method = radi.treutler
        grids.prune = gen_grid.nwchem_prune
        grids.radii_adjust = radi.treutler_atomic_radii_adjust
        grids.atomic_radii = radi.BRAGG_RADII


def calculate_molecule(config, directory, progress=print):
    try:
        import pyscf
        from pyscf import dft, gto, lib
        from pyscf.tools import cubegen
    except ImportError as exc:
        raise RuntimeError('PySCF가 없습니다. START_HERE_KO.md의 DFT 엔진 설치를 진행하세요.') from exc
    started = time.perf_counter()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    save_dft_config(config, directory / 'settings.json')
    check_cancelled(directory)
    lib.num_threads(config.threads)
    progress(f'{config.formula} · 전자 {config.electrons}개 · {config.resolved_method.upper()} / '
             f'{config.xc} / {config.basis_label} · SCF {config.scf_solver}')
    if config.resolved_ecp:
        progress(f'def2-ECP: {", ".join(config.resolved_ecp)} · 코어 {config.ecp_core_electrons}개 제외 / '
                 f'명시적 계산 전자 {config.explicit_electrons}개')
    try:
        mol = gto.M(atom=[(a.element, a.position) for a in config.atoms],
                    unit=config.unit, charge=config.charge, spin=config.spin,
                    basis=config.basis_by_element, ecp=config.resolved_ecp,
                    nucmod='gaussian' if config.resolved_nuclear_model == 'gaussian' else 0,
                    verbose=0, max_memory=config.memory_mb)
    except Exception as exc:
        raise ValueError(f'원자·기저함수·전하/스핀 설정을 구성할 수 없습니다: {exc}') from exc
    atomic_numbers = np.asarray([ATOMIC_NUMBERS[a.element] for a in config.atoms], dtype=int)
    effective_charges = mol.atom_charges()
    core_counts = np.asarray([mol.atom_nelec_core(i) for i in range(mol.natm)], dtype=int)
    expected_cores = np.asarray([DEF2_CORE_ELECTRONS[a.element] if a.element in config.resolved_ecp else 0
                                for a in config.atoms], dtype=int)
    if (mol.nelectron != config.explicit_electrons
            or not np.array_equal(core_counts, expected_cores)
            or not np.array_equal(atomic_numbers - core_counts, effective_charges)):
        raise RuntimeError('설정의 전자 수·ECP 코어와 엔진의 실제 구성이 일치하지 않습니다.')
    if mol.nao_nr() > 800:
        raise ValueError('기저함수 800개를 넘는 계산은 현재 데스크톱 엔진에서 지원하지 않습니다.')
    if max(mol.nelec) > mol.nao_nr():
        raise ValueError('선택한 기저함수에 전자를 모두 배치할 수 없습니다. 전하와 기저함수를 확인하세요.')
    dispersion = None
    if config.dispersion == 'd4':
        # Before SCF so a missing library fails fast. Atomic numbers, not ECP effective charges.
        from .dispersion import d4_energy
        dispersion = d4_energy(atomic_numbers, mol.atom_coords(), config.charge, config.xc)
        progress(f'D4 분산 보정 ({dispersion["library"]} {dispersion["library_version"]}, {config.xc} 매개변수): '
                 f'{dispersion["energy_hartree"]:.10f} Ha · SCF 후 총에너지에 더합니다.')
    if config.engine == 'native':
        from .native_scf import NativeKS
        mf = NativeKS(mol, config, directory)
    else:
        mf = dft.RKS(mol) if config.resolved_method == 'rks' else dft.UKS(mol)
        mf.init_guess = '1e' if config.initial_guess == 'core' else 'minao'
        mf.small_rho_cutoff = 0  # Fixed grids make native/reference comparisons reproducible.
    progress(f'계산 엔진: {"QuantumStarter 자체 SCF" if config.engine == "native" else "외부 PySCF SCF"}')
    if config.resolved_relativity == 'x2c':
        mf = mf.sfx2c1e()
    progress(f'해밀토니안: {config.resolved_relativity} · 핵 모형: {config.resolved_nuclear_model}')
    mf.xc = config.xc
    configure_grids(mf.grids, config)
    progress(f'적분 격자: {config.resolved_grid_scheme} · level {config.grid_level}')
    if config.scf_solver == 'newton':
        mf = mf.newton()
    mf.max_cycle = config.max_cycle
    mf.conv_tol = config.conv_tol
    mf.conv_tol_grad = float(np.sqrt(config.conv_tol))
    mf.chkfile = None
    history = []

    def iteration(env):
        cycle = env['imacro'] if config.scf_solver == 'newton' else env['cycle']
        row = {'cycle': int(cycle) + 1, 'solver': config.scf_solver,
               'energy_hartree': float(env['e_tot']),
               'gradient_norm': float(env['norm_gorb']) if env.get('norm_gorb') is not None else None,
               'density_change_norm': float(env['norm_ddm']) if env.get('norm_ddm') is not None else None,
               'delta_energy_hartree': float(env['delta_energy']) if env.get('delta_energy') is not None else None}
        checkpoint_dm = env.get('dm')
        if checkpoint_dm is None and env.get('mo_coeff') is not None and env.get('mo_occ') is not None:
            cs, occ = np.asarray(env['mo_coeff']), np.asarray(env['mo_occ'])
            checkpoint_dm = (cs * occ[..., None, :]) @ cs.swapaxes(-1, -2)
        if checkpoint_dm is not None:
            save_checkpoint(directory, config, checkpoint_dm, row['cycle'])
        check_cancelled(directory)
        # Newton reports its final macroiteration a second time after canonicalization.
        if history and history[-1]['cycle'] == row['cycle']:
            history[-1] = row
            return
        history.append(row)
        residual = f'{row["gradient_norm"]:.3e}' if row['gradient_norm'] is not None else '—'
        delta = f'{row["density_change_norm"]:.3e}' if row['density_change_norm'] is not None else '—'
        progress(f'SCF {row["cycle"]:3d} · E = {row["energy_hartree"]:.12f} Ha · 잔차 {residual} · ΔD {delta}')
        progress('SCF_JSON ' + json.dumps(row))

    mf.callback = iteration
    restart = directory / 'restart_input.npz'
    dm0 = load_checkpoint(restart, config, mol.intor_symmetric('int1e_ovlp')) if restart.exists() else None
    restart_engine = None
    if dm0 is not None:
        with np.load(restart, allow_pickle=False) as checkpoint:
            restart_engine = str(checkpoint['source_engine'])
        progress('호환성 검사를 통과한 밀도 체크포인트에서 재시작합니다. DIIS 이력은 새로 시작합니다.')
    energy = float(mf.kernel(dm0=dm0))
    if not np.isfinite(energy):
        raise RuntimeError('유한한 SCF 에너지를 얻지 못했습니다.')
    dm = np.asarray(mf.make_rdm1())
    save_checkpoint(directory, config, dm, len(history))
    check_cancelled(directory)
    coeff = np.asarray(mf.mo_coeff)
    energies = np.asarray(mf.mo_energy)
    occupations = np.asarray(mf.mo_occ)
    overlap = mol.intor_symmetric('int1e_ovlp')
    total_dm = dm.sum(axis=0) if dm.ndim == 3 else dm
    electron_trace = float(np.einsum('ij,ji->', total_dm, overlap))
    grid_density = mf._numint.get_rho(mol, total_dm, mf.grids, max_memory=config.memory_mb)
    grid_electrons = float(np.dot(grid_density, mf.grids.weights))
    spin_trace = float(np.einsum('ij,ji->', dm[0] - dm[1], overlap)) if dm.ndim == 3 else 0.0
    orbitals = []
    channels = [('restricted', energies, occupations, coeff)] if dm.ndim == 2 else [
        ('alpha', energies[0], occupations[0], coeff[0]),
        ('beta', energies[1], occupations[1], coeff[1])]
    warnings = []
    if abs(grid_electrons - electron_trace) > 0.01:
        warnings.append(f'DFT 적분 격자의 전자 수 오차 {grid_electrons - electron_trace:+.6f}개. '
                        'grid_level을 높여 수렴을 확인하세요. SCF 수렴과 격자 정확도는 별개입니다.')
    if config.resolved_ecp:
        warnings.append(f'def2-ECP 적용: 전체 전자 {config.electrons}개 중 코어 {config.ecp_core_electrons}개를 '
                        f'제외한 {config.explicit_electrons}개를 계산합니다. 저장 밀도에는 코어 전자밀도가 없으며, '
                        '전체전자 계산과 절대 총에너지를 직접 비교하지 마세요.')
    elif any(z > 36 for z in atomic_numbers) and config.resolved_relativity == 'none':
        warnings.append('무거운 원소의 전체전자 계산에 비상대론 근사를 선택했습니다. 상대론 효과 검증이 필요합니다.')
    if config.resolved_relativity == 'x2c':
        warnings.append('스칼라 X2C의 일전자 상대론 보정입니다. SOC와 이전자 picture-change 보정은 포함하지 않습니다. '
                        '밀도·오비탈 그림은 X2C 표현이며 물리적 사성분 밀도의 picture-change 보정은 하지 않았습니다.')
    if any(z >= 89 for z in atomic_numbers):
        warnings.append('악티늄족·초중원소: 기저 데이터와 계산 경로의 지원이며 바닥상태·SOC·정량 정확도를 보증하지 않습니다.')
    transition_metals = sorted({atom.element for atom in config.atoms} & TRANSITION_METALS)
    if transition_metals:
        warnings.append(f'전이금속 {", ".join(transition_metals)} 포함: 여러 2S와 초기 전자배치의 해가 '
                        '가능합니다. SCF 수렴만으로 바닥상태가 확인되지 않으며 스핀 상태 비교가 필요합니다.')
    if not mf.converged:
        warnings.append('SCF 미수렴: 마지막 반복의 진단값입니다. 성공한 계산으로 해석하지 마세요.')
        if config.scf_solver == 'diis' and config.engine == 'pyscf':
            warnings.append('SCF 반복 방식을 newton으로 바꾸거나 입력 스핀·기저·반복 제한을 점검하세요.')
    bases_used = set(config.basis_by_element.values())
    if 'sto-3g' in bases_used:
        warnings.append('STO-3G는 최소 기저입니다. 정량 비교 전 더 큰 기저로 수렴을 확인하세요.')
    if config.charge < 0 and not any('+' in basis for basis in bases_used):
        warnings.append('음이온은 확산함수가 중요합니다. 현재 선택한 기저에는 확산함수가 없습니다.')
    if len(bases_used) > 1:
        warnings.append(f'원소별 기저 혼합: {config.basis_label}. 다른 기저 조합의 총에너지와 직접 비교하지 말고, '
                        '반응·흡착에너지는 모든 항에서 원소별 기저가 같을 때만 계산하세요.')
        if config.resolved_relativity == 'x2c' and bases_used - {'dyall-v2z', 'dyall-v3z'}:
            warnings.append('X2C 계산에 비상대론용으로 수축된 기저가 섞여 있습니다. '
                            '해당 원소의 상대론 기저 수렴을 별도로 확인하세요.')
    if dispersion is not None:
        warnings.append(f'D4 분산 보정 {dispersion["energy_hartree"]:+.10f} Ha가 총에너지에 포함되었습니다. '
                        '밀도와 무관한 기하 보정이라 SCF 이력·밀도·오비탈에는 들어 있지 않습니다. '
                        '분산 보정이 다른 계산과 에너지를 섞지 마세요.')
        if config.xc == 'b3lyp':
            warnings.append('B3LYP는 PySCF/Libxc 정의(Libxc ID 402)입니다. D4 매개변수 피팅에 쓰인 '
                            'B3LYP 구현과 VWN 변형이 다를 수 있습니다.')
    if config.unused_basis_overrides:
        warnings.append('구조에 없는 원소의 기저 지정은 사용되지 않았습니다: '
                        + ', '.join(config.unused_basis_overrides))
    warnings.append('SCF 수렴은 정확도 보증이 아닙니다. 기저·격자·스핀 상태를 바꿔 확인하세요.')
    if abs(electron_trace - config.explicit_electrons) > 1e-5:
        raise RuntimeError('밀도행렬의 전자 수가 명시적으로 계산하는 전자 수와 일치하지 않습니다.')
    if abs(spin_trace - config.spin) > 1e-5:
        raise RuntimeError('밀도행렬의 스핀 전자 수 차이가 입력 2S와 일치하지 않습니다.')
    progress('전자밀도와 오비탈 데이터를 저장합니다.')
    cube_args = dict(nx=config.cube_points, ny=config.cube_points, nz=config.cube_points,
                     margin=config.cube_margin_bohr)
    density = cubegen.density(mol, str(directory / 'density.cube'), total_dm, **cube_args)
    cube = cubegen.Cube(mol, **cube_args)
    np.savez_compressed(directory / 'density_grid.npz', density=density,
                        origin_bohr=cube.boxorig, axes_bohr=cube.box,
                        atom_coords_bohr=mol.atom_coords(), atom_numbers=atomic_numbers,
                        effective_nuclear_charges=effective_charges, ecp_core_electrons_by_atom=core_counts)
    np.savez_compressed(directory / 'electronic_structure.npz', density_matrix=dm,
                        mo_coeff=coeff, mo_energy_hartree=energies, mo_occupation=occupations,
                        overlap=overlap, atom_coords_bohr=mol.atom_coords(),
                        atom_numbers=atomic_numbers, effective_nuclear_charges=effective_charges,
                        ecp_core_electrons_by_atom=core_counts, ao_labels=np.asarray(mol.ao_labels()))
    cube_files = ['density.cube']
    if dm.ndim == 3:
        cubegen.density(mol, str(directory / 'spin_density.cube'), dm[0] - dm[1], **cube_args)
        cube_files.append('spin_density.cube')
    for channel, es, occs, cs in channels:
        occupied = np.flatnonzero(occs > 1e-8)
        empty = np.flatnonzero(occs <= 1e-8)
        frontier = {}
        if len(occupied):
            frontier[int(occupied[-1])] = 'HOMO'
        if len(empty):
            frontier[int(empty[0])] = 'LUMO'
        for index, (e, occ) in enumerate(zip(es, occs)):
            cube_name = None
            label = frontier.get(index, '')
            if label:
                cube_name = f'{channel}_{label.lower()}.cube'
                cubegen.orbital(mol, str(directory / cube_name), cs[:, index], **cube_args)
                cube_files.append(cube_name)
            orbitals.append({'channel': channel, 'index': index + 1, 'label': label,
                             'occupation': float(occ), 'energy_hartree': float(e),
                             'energy_ev': float(e * HARTREE_TO_EV), 'cube_file': cube_name})
    if config.resolved_ecp or config.resolved_relativity == 'x2c':
        # Make the density model visible even when a CUBE is separated from its JSON.
        for name in cube_files:
            path = directory / name
            with path.open('r+', encoding='utf-8') as stream:
                lines = stream.readlines()
                lines[0] = ('ECP model; core electrons omitted. ' if config.resolved_ecp else
                            'Scalar X2C representation; no picture-change correction. ') + lines[0]
                stream.seek(0)
                stream.writelines(lines)
                stream.truncate()
    with (directory / 'orbitals.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(orbitals[0]))
        writer.writeheader()
        writer.writerows(orbitals)
    with (directory / 'scf.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=['cycle', 'solver', 'energy_hartree', 'gradient_norm', 'density_change_norm', 'delta_energy_hartree'])
        writer.writeheader()
        writer.writerows(history)
    spin_square = float(mf.spin_square()[0]) if dm.ndim == 3 else 0.0
    target_s2 = config.spin / 2 * (config.spin / 2 + 1)
    if abs(spin_square - target_s2) > 0.1:
        warnings.append('UKS 오비탈의 <S²>가 지정 스핀의 S(S+1)과 다릅니다. 스핀 오염을 점검하세요.')
    components = (mf.components if config.engine == 'native' else {
        'one_electron': float(np.einsum('ij,ji->', total_dm, mf.get_hcore())),
        'coulomb': float(mf.scf_summary['coul']), 'exchange_correlation': float(mf.scf_summary['exc']),
        'nuclear_repulsion': float(mol.energy_nuc())})
    scf_energy = energy
    if dispersion is not None:
        components = dict(components, dispersion=dispersion['energy_hartree'])
        energy = scf_energy + dispersion['energy_hartree']
    if abs(sum(components.values()) - energy) > 1e-7:
        raise RuntimeError('에너지 성분 합과 총에너지가 일치하지 않습니다.')
    xc_hybrid, xc_terms = dft.libxc.parse_xc(config.xc)
    if dm0 is not None:
        actual_guess = 'checkpoint'
    elif config.initial_guess != 'auto':
        actual_guess = config.initial_guess
    else:
        actual_guess = 'core' if config.engine == 'native' else 'minao'
    summary = {
        'schema_version': 4, 'program': 'QuantumStarter', 'version': '0.8.2',
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'engine': 'quantumstarter-native' if config.engine == 'native' else 'pyscf-molecular', 'system': 'molecule',
        'provenance': {'scf_driver': 'QuantumStarter' if config.engine == 'native' else 'PySCF',
                       'integrals': 'PySCF/libcint', 'xc_library': 'Libxc',
                       'xc_functional_ids_and_weights': [[int(i), float(w)] for i, w in xc_terms],
                       'exact_exchange_fraction': float(dft.libxc.hybrid_coeff(config.xc)),
                       'restarted': dm0 is not None, 'restart_source_engine': restart_engine},
        'initial_guess_details': {'requested': config.initial_guess,
                                  'resolved': actual_guess,
                                  'perturbation_seed': 1729 if config.engine == 'native' and config.initial_guess == 'perturbed' and dm0 is None else None},
        'energy_components_hartree': components, 'energy_component_sum_error_hartree': sum(components.values()) - energy,
        'scf_criteria': ({'energy_hartree': config.conv_tol, 'density_frobenius': config.native_density_tol,
                          'orthogonal_commutator_frobenius': config.native_residual_tol,
                          'final_residual': mf.final_residual} if config.engine == 'native' else
                         {'energy_hartree': config.conv_tol, 'orbital_gradient': mf.conv_tol_grad}),
        'checkpoint_file': 'checkpoint.npz',
        'status': 'converged' if mf.converged else 'not_converged',
        'converged': bool(mf.converged), 'config': config.to_dict(), 'formula': config.formula,
        'method': config.resolved_method.upper(), 'xc': config.xc, 'basis': config.basis,
        'basis_by_element': config.basis_by_element,
        'scf_solver': config.scf_solver,
        'relativity': config.resolved_relativity, 'nuclear_model': config.resolved_nuclear_model,
        'integration_grid': {'scheme': config.resolved_grid_scheme, 'level': config.grid_level,
                             'atom_grid_overrides': {key: [int(n) for n in value]
                                                     for key, value in mf.grids.atom_grid.items()},
                             'radial_method': mf.grids.radi_method.__name__,
                             'pruning': mf.grids.prune.__name__ if mf.grids.prune else None,
                             'radii_adjustment': mf.grids.radii_adjust.__name__ if mf.grids.radii_adjust else None,
                             'points': len(mf.grids.weights), 'integrated_electrons': grid_electrons,
                             'electron_count_error': grid_electrons - electron_trace},
        'nuclear_gaussian_zeta_per_bohr2': [float(mol._env[mol._atm[i, gto.PTR_ZETA]]) for i in range(mol.natm)],
        'scf_iteration_kind': 'macroiteration' if config.scf_solver == 'newton' else 'iteration',
        'total_energy_hartree': energy, 'total_energy_ev': energy * HARTREE_TO_EV,
        'scf_energy_hartree': scf_energy, 'dispersion': dispersion,
        'energy_definition': ('SCF + D4 dispersion; scf_history holds SCF energies only'
                              if dispersion is not None else 'SCF'),
        'nuclear_repulsion_hartree': float(mol.energy_nuc()),
        'electrons': config.electrons, 'electron_count_from_density': electron_trace,
        'explicit_electrons': config.explicit_electrons, 'ecp_core_electrons': config.ecp_core_electrons,
        'ecp': {'requested': config.ecp, 'resolved': config.resolved_ecp,
                'core_electrons_by_atom': core_counts.tolist(),
                'effective_nuclear_charges': effective_charges.tolist()},
        'atomic_numbers': atomic_numbers.tolist(),
        'spin_from_density': spin_trace,
        'alpha_electrons': int(mol.nelec[0]), 'beta_electrons': int(mol.nelec[1]),
        'basis_functions': int(mol.nao_nr()), 'spin_square': spin_square,
        'target_spin_square': target_s2, 'spin_square_deviation': spin_square - target_s2,
        'scf_cycles': len(history), 'scf_history': history, 'orbitals': orbitals,
        'density': {'file': 'density.cube', 'grid_file': 'density_grid.npz',
                    'representation': ('ecp_pseudo_density' if config.resolved_ecp else
                                       'x2c_picture_density' if config.resolved_relativity == 'x2c' else 'all_electron_density'),
                    'picture_change_corrected': False,
                    'includes_core_electrons': not bool(config.resolved_ecp),
                    'electron_count': config.explicit_electrons,
                    'unit': 'electron/bohr^3', 'points_per_axis': config.cube_points,
                    'margin_bohr': config.cube_margin_bohr},
        'orbital_data_file': 'electronic_structure.npz', 'cube_files': cube_files,
        'units': {'energy': 'hartree', 'coordinates_in_cubes': 'bohr',
                  'density': 'electron/bohr^3', 'orbital_amplitude': 'bohr^(-3/2)'},
        'assumptions': ['isolated finite molecule', 'fixed nuclei, single point energy',
                        'scalar def2-ECP KS-DFT; no spin-orbit coupling' if config.resolved_ecp
                        else 'scalar X2C all-electron KS-DFT; no SOC' if config.resolved_relativity == 'x2c'
                        else 'nonrelativistic all-electron KS-DFT', 'no periodic boundary conditions']
                       + (['DFT-D4 dispersion (EEQ charges, BJ damping, ATM three-body) added after SCF']
                          if dispersion is not None else []),
        'warnings': warnings, 'elapsed_seconds': time.perf_counter() - started,
        'environment': {'python': platform.python_version(), 'pyscf': pyscf.__version__,
                        'numpy': np.__version__, 'libxc': dft.libxc.__version__},
    }
    # results.json is the commit marker: never publish success before all raw files exist.
    temporary = directory / 'results.json.tmp'
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(directory / 'results.json')
    progress('SCF 수렴 · 결과 저장 완료' if mf.converged else 'SCF 미수렴 · 진단 데이터 저장 완료')
    return summary
