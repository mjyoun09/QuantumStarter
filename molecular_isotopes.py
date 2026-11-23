"""Compare harmonic vibrations after mass substitution on a saved molecular Hessian."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from engines.frequencies import (BOHR_ANGSTROM, HARTREE_TO_KJ_MOL, STATIONARY_GRADIENT_LIMIT,
                                 harmonic_analysis, zero_point_energy_hartree)
from engines.gas_thermo import write_csv
from molecular import Atom
from workflow import VERSION, new_run_dir

ROOT = Path(__file__).resolve().parent
WATER_EXAMPLE = ROOT / 'examples' / 'water_frequency_reference.json'
# NIST, Atomic Weights and Isotopic Compositions for Hydrogen (neutral ²H, u).
DEUTERIUM_MASS_AMU = 2.01410177812


def substituted_masses(elements, reference, deuterate, assignments):
    if bool(deuterate) == bool(assignments):
        raise ValueError('수소 전체 치환 또는 원자별 질량 지정 중 하나를 선택하세요.')
    masses = reference.copy()
    if deuterate:
        hydrogen = np.asarray(elements) == 'H'
        if not hydrogen.any():
            raise ValueError('중수소로 바꿀 수소 원자가 없습니다.')
        masses[hydrogen] = DEUTERIUM_MASS_AMU
        return masses

    seen = set()
    for token in assignments:
        try:
            number, value = token.split('=')
            atom, mass = int(number), float(value)
        except ValueError:
            raise ValueError(f'원자별 질량은 2=2.0141처럼 원자 순번=질량(u)으로 쓰세요: {token}') from None
        if not 1 <= atom <= len(elements) or atom in seen:
            raise ValueError(f'원자 순번은 1–{len(elements)} 사이에서 한 번씩 지정하세요: {atom}')
        if not np.isfinite(mass) or mass <= 0:
            raise ValueError(f'{atom}번 원자의 질량은 유한한 양수여야 합니다.')
        masses[atom - 1] = mass
        seen.add(atom)
    return masses


def compare_isotopes(summary, *, deuterate=False, assignments=()):
    if (not isinstance(summary, dict) or summary.get('system', 'molecule') != 'molecule'
            or summary.get('calculation') != 'frequencies' or summary.get('converged') is not True):
        raise ValueError('수렴한 분자 진동수 계산의 results.json을 선택하세요.')
    try:
        config = summary['config']
        atoms = [Atom(**item) for item in config['atoms']]
        unit = config['unit']
        vibrations = summary['frequencies']
        reference = np.asarray(vibrations['masses_amu'], dtype=float)
        hessian = np.asarray(vibrations['hessian_hartree_per_bohr2'], dtype=float)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f'결과 파일의 좌표·질량·헤시안을 읽을 수 없습니다 ({exc}).') from exc
    n = len(atoms)
    if n < 2 or unit not in ('angstrom', 'bohr'):
        raise ValueError('좌표 단위(angstrom 또는 bohr)와 원자 두 개 이상의 구조가 필요합니다.')
    if reference.shape != (n,) or not np.isfinite(reference).all() or np.any(reference <= 0):
        raise ValueError('기준 질량은 원자마다 하나씩, 유한한 양수여야 합니다.')
    if (hessian.shape != (3 * n, 3 * n) or not np.isfinite(hessian).all()
            or not np.allclose(hessian, hessian.T, atol=1e-10, rtol=1e-10)):
        raise ValueError('결과에 유한한 대칭 3N × 3N 헤시안(Ha/bohr²)이 필요합니다.')

    elements = [atom.element for atom in atoms]
    masses = substituted_masses(elements, reference, deuterate, assignments)
    coords = np.asarray([atom.position for atom in atoms])
    if unit == 'angstrom':
        coords = coords / BOHR_ANGSTROM
    original = harmonic_analysis(hessian, coords, reference)
    changed = harmonic_analysis(hessian, coords, masses)
    if original['external_mode_count'] != changed['external_mode_count']:
        raise ValueError('치환 전후 병진·회전 모드 개수가 달라졌습니다. 구조의 선형성과 질량을 확인하세요.')

    # Compare displacements in the same Cartesian frame and reference mass metric.
    # Eigenvector signs are arbitrary. Maximize the total absolute overlap one to one.
    def unit_modes(analysis):
        modes = analysis['normal_modes'].reshape(-1, 3 * n) * np.repeat(np.sqrt(reference), 3)
        return modes / np.linalg.norm(modes, axis=1, keepdims=True)

    overlap = np.clip(np.abs(unit_modes(original) @ unit_modes(changed).T), 0, 1)
    before, after = linear_sum_assignment(-overlap)
    freqs = original['frequencies_cm1']
    new_freqs = changed['frequencies_cm1']
    rows = [dict(reference_mode=int(i + 1), substituted_mode=int(j + 1), reference_cm1=float(freqs[i]),
                 substituted_cm1=float(new_freqs[j]), shift_cm1=float(new_freqs[j] - freqs[i]),
                 absolute_overlap=float(overlap[i, j])) for i, j in zip(before, after)]
    zpe = zero_point_energy_hartree(freqs) * HARTREE_TO_KJ_MOL
    new_zpe = zero_point_energy_hartree(new_freqs) * HARTREE_TO_KJ_MOL
    gradient = vibrations.get('reference_max_gradient_hartree_per_bohr')
    warnings = []
    if gradient is not None and (not isinstance(gradient, (float, int)) or not np.isfinite(gradient) or gradient < 0):
        raise ValueError('기준 구조의 최대 기울기는 유한한 음이 아닌 수여야 합니다.')
    if gradient is not None and gradient > STATIONARY_GRADIENT_LIMIT:
        warnings.append(f'기준 구조의 최대 기울기가 {gradient:.2g} Ha/bohr입니다. 평형 구조인지 확인하세요.')
    if np.any(freqs < 0) or np.any(new_freqs < 0):
        warnings.append('허수 모드는 음수 파수로 표시했습니다. ZPE에는 양의 파수만 포함했습니다.')
    return {
        'version': VERSION, 'model': 'Born-Oppenheimer harmonic vibrations at the source geometry and Hessian',
        'matching': 'maximum total absolute Cartesian mode overlap in the reference mass metric',
        'formula': summary.get('formula'), 'engine': summary.get('engine'),
        'reference_max_gradient_hartree_per_bohr': gradient,
        'atoms': [dict(atom=i, element=e, reference_mass_amu=float(a), substituted_mass_amu=float(b))
                  for i, (e, a, b) in enumerate(zip(elements, reference, masses), 1)],
        'modes': rows, 'reference_zpe_kj_mol': zpe, 'substituted_zpe_kj_mol': new_zpe,
        'delta_zpe_kj_mol': new_zpe - zpe,
        'warnings': warnings, 'source_warnings': summary.get('warnings', []),
    }


def run_isotopes(source, output_root=None, *, deuterate=False, assignments=()):
    source = Path(source).resolve()
    if source.is_dir():
        source = source / 'results.json'
    raw = source.read_bytes()
    report = compare_isotopes(json.loads(raw), deuterate=deuterate, assignments=assignments)
    report.update(source_file=str(source), source_sha256=hashlib.sha256(raw).hexdigest())
    directory = new_run_dir(output_root or ROOT / 'outputs' / 'isotopes')
    write_csv(directory / 'modes.csv', report['modes'])
    (directory / 'isotopes.json').write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    return report, directory


def main(argv=None):
    parser = argparse.ArgumentParser(prog='main.py isotopes',
                                     description='저장된 분자 헤시안으로 동위원소 진동수와 영점에너지를 비교합니다.')
    parser.add_argument('source', type=Path, help='진동수 results.json 또는 해당 실행 폴더')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--deuterate', action='store_true', help='모든 수소를 중수소(²H)로 치환')
    group.add_argument('--mass', nargs='+', metavar='INDEX=AMU', help='입력 순서(1부터)와 질량(u), 예: 2=2.0141 3=2.0141')
    parser.add_argument('--output', type=Path, help='새 비교 결과 폴더를 만들 상위 경로')
    args = parser.parse_args(argv)
    try:
        report, directory = run_isotopes(args.source, args.output, deuterate=args.deuterate, assignments=args.mass)
    except (ValueError, OSError) as exc:
        parser.exit(1, f'동위원소 비교 오류: {exc}\n')
    for row in report['modes']:
        print(f"{row['reference_mode']} → {row['substituted_mode']}    {row['reference_cm1']:.2f} → "
              f"{row['substituted_cm1']:.2f} cm⁻¹    Δ {row['shift_cm1']:+.2f} cm⁻¹")
    print(f"ZPE {report['reference_zpe_kj_mol']:.3f} → {report['substituted_zpe_kj_mol']:.3f} kJ/mol "
          f"(Δ {report['delta_zpe_kj_mol']:+.3f})\n{directory}")
    for warning in report['warnings']:
        print(warning)
    return 0
