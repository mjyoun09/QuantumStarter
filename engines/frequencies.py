"""Harmonic vibrational analysis from central finite differences of analytic gradients.

H[:, j] = (g(x + h e_j) - g(x - h e_j)) / 2h, symmetrized. Masses are the most
abundant isotopes (exact isotopic masses). Translations and rotations (5 for
linear molecules) are removed exactly by diagonalizing the mass-weighted Hessian
in the orthogonal complement of their span, so 3N-6 (3N-5) frequencies remain.
Every displaced calculation is a full recorded run under displacements/.
"""
from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, Boltzmann, Planck, physical_constants, speed_of_light

from molecular import Atom, save_dft_config
from scf_control import check_cancelled
from .pyscf_engine import BOHR_ANGSTROM, HARTREE_TO_EV, calculate_molecule

HARTREE_J = physical_constants['Hartree energy'][0]
BOHR_M = physical_constants['Bohr radius'][0]
AMU_KG = physical_constants['atomic mass constant'][0]
# sqrt(Hartree / (bohr^2 amu)) in rad/s -> wavenumber in cm^-1
AU_TO_CM1 = np.sqrt(HARTREE_J / (BOHR_M**2 * AMU_KG)) / (2 * np.pi * speed_of_light * 100)
CM1_TO_HARTREE = Planck * speed_of_light * 100 / HARTREE_J
HARTREE_TO_KJ_MOL = HARTREE_J * Avogadro / 1000
SIGNIFICANT_IMAGINARY_CM1 = 20.0


def isotope_masses(elements):
    from pyscf.data.elements import COMMON_ISOTOPE_MASSES, ELEMENTS
    masses = [float(COMMON_ISOTOPE_MASSES[ELEMENTS.index(e)]) for e in elements]
    if min(masses) <= 0:
        raise ValueError('일부 원소의 동위원소 질량 데이터가 없습니다. 진동수를 계산할 수 없습니다.')
    return np.asarray(masses)


def harmonic_analysis(hessian, coords_bohr, masses):
    """Frequencies [cm^-1, imaginary as negative], Cartesian normal modes and diagnostics."""
    natm = len(masses)
    root = np.repeat(np.sqrt(masses), 3)
    weighted = hessian / np.outer(root, root)
    centred = coords_bohr - masses @ coords_bohr / masses.sum()
    rigid = []
    for axis in range(3):
        translation = np.zeros((natm, 3))
        translation[:, axis] = 1
        rotation = np.cross(np.eye(3)[axis], centred)
        rigid += [translation.ravel() * root, rotation.ravel() * root]
    u, s, _ = np.linalg.svd(np.array(rigid).T, full_matrices=True)
    rank = int(np.sum(s > 1e-6 * s.max()))
    external, internal = u[:, :rank], u[:, rank:]
    values, vectors = np.linalg.eigh(internal.T @ weighted @ internal)
    frequencies = np.sign(values) * np.sqrt(np.abs(values)) * AU_TO_CM1
    cartesian = (internal @ vectors) / root[:, None]
    modes = (cartesian / np.linalg.norm(cartesian, axis=0)).T.reshape(-1, natm, 3)
    residual = np.linalg.eigvalsh(external.T @ weighted @ external)
    return {'frequencies_cm1': frequencies, 'normal_modes': modes, 'linear': rank == 5,
            'external_mode_count': rank,
            'translation_rotation_residual_cm1': np.sign(residual) * np.sqrt(np.abs(residual)) * AU_TO_CM1}


def vibrational_thermo(frequencies_cm1, temperature):
    """Harmonic-oscillator ZPE and thermal vibrational terms for real modes; Hartree per molecule."""
    real = np.asarray([f for f in frequencies_cm1 if f > 0])
    zpe = 0.5 * float(np.sum(real)) * CM1_TO_HARTREE
    x = real * CM1_TO_HARTREE * HARTREE_J / (Boltzmann * temperature)  # h c nu / k T
    kt = Boltzmann * temperature / HARTREE_J
    # Written with q = exp(-x), which underflows to 0 instead of overflowing at low T.
    q = np.exp(-x)
    occupation = q / -np.expm1(-x)  # 1 / (e^x - 1)
    thermal = float(np.sum(real * CM1_TO_HARTREE * occupation))
    entropy = float(np.sum(x * occupation - np.log1p(-q))) * Boltzmann / HARTREE_J
    heat_capacity = float(np.sum(x**2 * occupation * (1 + occupation))) * Boltzmann / HARTREE_J
    return {'temperature_k': temperature, 'zpe_hartree': zpe, 'zpe_kj_mol': zpe * HARTREE_TO_KJ_MOL,
            'zpe_ev': zpe * HARTREE_TO_EV, 'u_vib_hartree': zpe + thermal,
            's_vib_hartree_per_k': entropy, 'cv_vib_hartree_per_k': heat_capacity,
            'minus_t_s_vib_hartree': -temperature * entropy, 'kt_hartree': kt,
            'model': 'harmonic oscillator, real modes only; translation/rotation/electronic terms not included'}


def symmetry_number(elements, coords_bohr):
    """Rotational symmetry number from PySCF's point-group detection (external helper)."""
    from pyscf import gto
    from pyscf.hessian.thermo import rotational_symmetry_number
    mol = gto.Mole()
    mol._atom = [(e, tuple(c)) for e, c in zip(elements, coords_bohr)]
    return rotational_symmetry_number(mol)


def ideal_gas_thermo(electronic_energy, frequencies_cm1, masses, coords_bohr, multiplicity, sigma,
                     temperature, pressure):
    """Rigid-rotor / harmonic-oscillator ideal gas; energies in Hartree, entropies in Hartree/K."""
    kb = Boltzmann / HARTREE_J
    kt = kb * temperature
    vib = vibrational_thermo(frequencies_cm1, temperature)
    total_mass = masses.sum() * AMU_KG
    q_trans = (2 * np.pi * total_mass * Boltzmann * temperature / Planck**2)**1.5 * Boltzmann * temperature / pressure
    centred = coords_bohr - masses @ coords_bohr / masses.sum()
    inertia = np.einsum('i,ij,ik->jk', masses, centred, centred)
    moments = np.sort(np.linalg.eigvalsh(np.trace(inertia) * np.eye(3) - inertia))  # amu bohr^2
    def theta(moment):  # Rotational temperature h^2 / (8 pi^2 I k_B) of a nonzero moment.
        return Planck**2 / (8 * np.pi**2 * moment * AMU_KG * BOHR_M**2 * Boltzmann)
    if len(masses) == 1:
        rotor, s_rot, u_rot = 'atom', 0.0, 0.0
    elif moments[0] < 1e-6 * moments[2]:
        rotor = 'linear'
        s_rot, u_rot = kb * (1 + np.log(temperature / (sigma * theta(moments[2])))), kt
    else:
        rotor = 'nonlinear'
        q_rot = np.sqrt(np.pi) / sigma * np.sqrt(temperature**3 / np.prod([theta(m) for m in moments]))
        s_rot, u_rot = kb * (1.5 + np.log(q_rot)), 1.5 * kt
    terms = {
        'electronic': {'energy': electronic_energy, 'entropy': kb * np.log(multiplicity)},
        'translational': {'energy': 1.5 * kt, 'entropy': kb * (2.5 + np.log(q_trans))},
        'rotational': {'energy': u_rot, 'entropy': float(s_rot)},
        'vibrational': {'energy': vib['u_vib_hartree'], 'entropy': vib['s_vib_hartree_per_k']}}
    internal = sum(t['energy'] for t in terms.values())
    entropy = sum(t['entropy'] for t in terms.values())
    enthalpy = internal + kt  # H = U + pV = U + kT for an ideal gas
    return {'model': 'ideal gas, rigid rotor, harmonic oscillator (real modes); no quasi-RRHO, no scaling',
            'temperature_k': temperature, 'pressure_pa': pressure, 'multiplicity': multiplicity,
            'symmetry_number': sigma, 'rotor': rotor, 'moments_of_inertia_amu_bohr2': moments.tolist(),
            'terms_hartree': terms, 'zpe_hartree': vib['zpe_hartree'],
            'zpe_corrected_energy_hartree': electronic_energy + vib['zpe_hartree'],
            'internal_energy_hartree': internal, 'enthalpy_hartree': enthalpy,
            'entropy_hartree_per_k': entropy, 'gibbs_free_energy_hartree': enthalpy - temperature * entropy,
            'thermal_correction_to_g_hartree': enthalpy - temperature * entropy - electronic_energy}


def frequency_molecule(config, directory, progress=print):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / 'restart_input.npz').exists():
        raise ValueError('진동수 계산은 밀도 체크포인트 재시작을 지원하지 않습니다. 새 실행으로 시작하세요.')
    save_dft_config(config, directory / 'settings.json')
    scale = 1.0 if config.unit == 'angstrom' else BOHR_ANGSTROM
    atoms = tuple(Atom(a.element, tuple(x * scale for x in a.position)) for a in config.atoms)
    base = replace(config, atoms=atoms, unit='angstrom', calculation='gradient', cube_points=12)
    masses = isotope_masses([a.element for a in atoms])
    coords = np.array([a.position for a in atoms]) / BOHR_ANGSTROM
    step = config.freq_step_bohr
    root = directory / 'displacements'
    count = 6 * len(atoms)

    def quiet(message):
        check_cancelled(directory)
        if message.startswith(('계산 오류', 'SCF 미수렴')):
            progress(message)

    def gradient_at(positions, name, guess):
        moved = tuple(Atom(a.element, tuple(p * BOHR_ANGSTROM)) for a, p in zip(atoms, positions))
        summary = calculate_molecule(replace(base, atoms=moved), root / name, quiet, guess)
        if not summary['converged'] or summary['gradient'] is None:
            raise RuntimeError(f'{name} 변위 계산의 SCF가 수렴하지 않아 진동수 계산을 중단합니다. {root / name}을 확인하세요.')
        with np.load(root / name / 'electronic_structure.npz', allow_pickle=False) as data:
            density = np.asarray(data['density_matrix'])
        return np.asarray(summary['gradient']['gradient_hartree_per_bohr']).ravel(), density, summary

    progress(f'진동수: 해석적 힘의 중심 유한차분 · 변위 {step} bohr · 변위 계산 {count}회 + 기준 1회')
    reference_gradient, reference_density, _ = gradient_at(coords, 'reference', None)
    hessian = np.zeros((coords.size, coords.size))
    for index in range(coords.size):
        atom, axis = divmod(index, 3)
        columns = []
        for sign, label in ((1, 'plus'), (-1, 'minus')):
            check_cancelled(directory)
            positions = coords.copy()
            positions[atom, axis] += sign * step
            columns.append(gradient_at(positions, f'atom{atom + 1:03d}_{"xyz"[axis]}_{label}', reference_density)[0])
        hessian[:, index] = (columns[0] - columns[1]) / (2 * step)
        progress(f'진동수 변위 {2 * index + 2}/{count} · 원자 {atom + 1} {"xyz"[axis]}')
    asymmetry = float(np.abs(hessian - hessian.T).max())
    hessian = (hessian + hessian.T) / 2
    analysis = harmonic_analysis(hessian, coords, masses)
    frequencies = analysis['frequencies_cm1']
    thermo = vibrational_thermo(frequencies, config.thermo_temperature_k)
    progress('기준 구조를 요청한 출력 설정으로 다시 계산합니다.')
    summary = calculate_molecule(replace(base, cube_points=config.cube_points), directory, progress, reference_density)
    save_dft_config(config, directory / 'settings.json')
    final_gradient = np.asarray(summary['gradient']['gradient_hartree_per_bohr'])
    sigma = config.symmetry_number or symmetry_number([a.element for a in atoms], coords)
    gas = ideal_gas_thermo(summary['total_energy_hartree'], frequencies, masses, coords, config.spin + 1, sigma,
                           config.thermo_temperature_k, config.thermo_pressure_pa)
    gas['symmetry_number_source'] = 'user' if config.symmetry_number else 'pyscf point-group detection'
    low = [float(f) for f in frequencies if 0 < f < 100]
    imaginary = [float(f) for f in frequencies if f < 0]
    significant = [f for f in imaginary if -f > SIGNIFICANT_IMAGINARY_CM1]
    summary['frequencies'] = {
        'method': 'central finite differences of analytic gradients', 'step_bohr': step,
        'gradient_calculations': count + 1, 'directory': 'displacements', 'grid_response': False,
        'masses_amu': masses.tolist(), 'mass_model': 'most abundant isotope (pyscf COMMON_ISOTOPE_MASSES)',
        'linear': analysis['linear'], 'vibrational_mode_count': len(frequencies),
        'frequencies_cm1': frequencies.tolist(), 'imaginary_frequencies_cm1': imaginary,
        'significant_imaginary_threshold_cm1': SIGNIFICANT_IMAGINARY_CM1,
        'stationary_point_character': ('minimum' if not significant else
                                       'transition_state_candidate' if len(significant) == 1 else 'higher_order_saddle'),
        'normal_modes_cartesian': analysis['normal_modes'].tolist(),
        'translation_rotation_residual_cm1': analysis['translation_rotation_residual_cm1'].tolist(),
        'hessian_hartree_per_bohr2': hessian.tolist(), 'hessian_asymmetry_max': asymmetry,
        'reference_max_gradient_hartree_per_bohr': float(np.abs(final_gradient).max()),
        'reference_gradient_recompute_difference': float(np.abs(final_gradient.ravel() - reference_gradient).max()),
        'thermochemistry': thermo, 'ideal_gas_thermochemistry': gas,
        'completed_at_utc': datetime.now(timezone.utc).isoformat()}
    # Top-level fields that the reaction/adsorption analysis can select.
    summary['zpe_corrected_energy_hartree'] = gas['zpe_corrected_energy_hartree']
    summary['enthalpy_hartree'] = gas['enthalpy_hartree']
    summary['gibbs_free_energy_hartree'] = gas['gibbs_free_energy_hartree']
    summary['calculation'] = 'frequencies'
    summary['config'] = config.to_dict()
    summary['assumptions'] = ['isolated finite molecule', 'harmonic vibrations at the input geometry',
                              *summary['assumptions'][2:]]
    warnings = []
    if np.abs(final_gradient).max() > 4.5e-4:
        warnings.append(f'입력 구조의 최대 기울기 {np.abs(final_gradient).max():.1e} Ha/bohr: 정지점이 아닙니다. '
                        '조화 진동수와 ZPE는 최적화된 구조에서만 물리적으로 해석하세요.')
    if significant:
        warnings.append(f'허수 진동 {len(significant)}개(|ν| > {SIGNIFICANT_IMAGINARY_CM1:g} cm⁻¹): 최소점이 아닙니다.')
    elif imaginary:
        warnings.append(f'작은 허수 진동 {len(imaginary)}개(≤ {SIGNIFICANT_IMAGINARY_CM1:g} cm⁻¹): 수치 잡음일 수 있습니다. '
                        '격자 수준·SCF/최적화 수렴 기준을 높여 확인하세요.')
    if low:
        warnings.append(f'100 cm⁻¹ 미만 저주파 모드 {len(low)}개: 조화 진동자 엔트로피가 과대평가될 수 있습니다'
                        '(quasi-RRHO 보정은 적용하지 않음). G 해석에 주의하세요.')
    warnings.append('조화 근사 진동수입니다. 비조화성·보정 인자 없음, ZPE/열역학에서 허수 모드 제외. '
                    f'G는 이상기체·강체회전자·조화진동자 근사(T={config.thermo_temperature_k:g} K, '
                    f'p={config.thermo_pressure_pa:g} Pa, 대칭수 {sigma:g}, 다중도 {config.spin + 1}, '
                    '들뜬 전자 상태 무시)입니다.')
    summary['warnings'] = warnings + summary['warnings']
    temporary = directory / 'results.json.tmp'
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(directory / 'results.json')
    lowest = ', '.join(f'{f:.1f}' for f in frequencies[:3])
    progress(f'진동수 {len(frequencies)}개 · 최저 {lowest} cm⁻¹ · 허수 {len(imaginary)}개 · '
             f'ZPE {thermo["zpe_kj_mol"]:.3f} kJ/mol · 저장 완료')
    return summary
