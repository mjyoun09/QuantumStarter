"""Auditable reaction/adsorption energies with composition and method checks."""
import argparse
from collections import Counter
from datetime import datetime
import hashlib
import json
import math
import sys
from pathlib import Path
from tempfile import mkdtemp
from scipy.constants import physical_constants
from molecular import DFTConfig
from periodic import config_from_dict

HARTREE_EV = physical_constants['Hartree energy in eV'][0]
ENERGY_FIELDS = ('total_energy_hartree', 'free_energy_hartree', 'extrapolated_energy_hartree')


def read_result(path):
    path = Path(path)
    if path.is_dir():
        path = path / 'results.json'
    raw = path.read_bytes()
    result = json.loads(raw.decode('utf-8-sig'))
    config_from_dict(result['config'])
    return result, {'path': str(path.resolve()), 'sha256': hashlib.sha256(raw).hexdigest()}


def compatibility(result):
    cfg = config_from_dict(result['config'])
    engine = result.get('engine')
    if engine not in ('quantumstarter-native', 'pyscf-molecular', 'pyscf-pbc'):
        raise ValueError('검증되지 않은 엔진/ML/수소 검증기 에너지는 이 DFT 반응 계산에 섞을 수 없습니다.')
    common = dict(engine=engine, system=cfg.system, xc=cfg.xc,
                  grid_level=cfg.grid_level, conv_tol=cfg.conv_tol)
    if isinstance(cfg, DFTConfig):
        # Molecules in one reaction contain different elements, so the basis is
        # compared per element (with its ECP) in reaction_energy, not as one name.
        ecp = cfg.resolved_ecp
        common.update(relativity=cfg.resolved_relativity, nuclear_model=cfg.resolved_nuclear_model,
                      grid_scheme=cfg.resolved_grid_scheme,
                      libxc_version=result.get('environment', {}).get('libxc'),
                      xc_definition=result.get('provenance', {}).get('xc_functional_ids_and_weights'),
                      basis_by_element={e: b + ('+def2-ecp' if e in ecp else '')
                                        for e, b in cfg.basis_by_element.items()})
    else:
        common.update(basis=cfg.basis, pseudo=cfg.pseudo, lattice=cfg.lattice, boundary=cfg.boundary,
                      kmesh=cfg.kmesh, smearing_ev=cfg.smearing_ev, precision=cfg.precision)
        if result.get('energy_normalization') != 'per_input_cell':
            raise ValueError('주기계 에너지는 입력 셀당 값이어야 합니다.')
    return common


def reaction_energy(terms, energy_field='total_energy_hartree'):
    """terms: sequence of (signed stoichiometric coefficient, result path).

    Products are positive, reactants negative. Geometry and spin may differ;
    composition and charge balance, Hamiltonian, grids and convergence may not.
    """
    if energy_field not in ENERGY_FIELDS or len(terms) < 2:
        raise ValueError('지원되는 에너지 종류와 두 개 이상의 반응 항이 필요합니다.')
    balance, charge, energy = Counter(), 0., 0.
    rows, reference, element_bases = [], None, {}
    for coefficient, path in terms:
        if isinstance(coefficient, bool) or not math.isfinite(coefficient) or coefficient == 0:
            raise ValueError('반응 계수는 0이 아닌 유한한 수여야 합니다.')
        result, source = read_result(path)
        if result.get('status') != 'converged' or result.get('converged') is not True:
            raise ValueError(f'수렴하지 않은 계산은 에너지 분석에 사용할 수 없습니다: {path}')
        key = compatibility(result)
        for element, basis in key.pop('basis_by_element', {}).items():
            if element_bases.setdefault(element, basis) != basis:
                raise ValueError(f'계산 조건 불일치: 원소 {element}의 기저/ECP가 항마다 다릅니다 '
                                 f'({element_bases[element]} ≠ {basis}). 같은 원소에는 같은 기저를 쓰세요.')
        if reference is not None:
            differences = sorted(k for k in reference.keys() | key.keys() if reference.get(k) != key.get(k))
            if differences:
                raise ValueError('계산 조건 불일치: ' + ', '.join(differences) + '. 같은 조건으로 계산하세요.')
        reference = key
        cfg = config_from_dict(result['config'])
        for atom in cfg.atoms:
            balance[atom.element] += coefficient
        charge += coefficient * cfg.charge
        value = result.get(energy_field)
        if not isinstance(value, (float, int)) or isinstance(value, bool) or not math.isfinite(value):
            raise ValueError(f'{path}: 선택한 에너지 값이 없거나 유한하지 않습니다.')
        energy += coefficient * value
        rows.append(dict(coefficient=coefficient, formula=cfg.formula, charge=cfg.charge, spin=cfg.spin,
                         engine=result['engine'], energy_hartree=value, source=source))
    imbalance = {e: v for e, v in balance.items() if abs(v) > 1e-8}
    if imbalance or abs(charge) > 1e-8:
        raise ValueError(f'원소/전하 수지가 맞지 않습니다: 원소 {imbalance}, 전하 {charge:+g}')
    if element_bases:
        reference = dict(reference, basis_by_element=dict(sorted(element_bases.items())))
    return {'schema_version': 1, 'analysis': 'reaction_energy', 'energy_field': energy_field,
            'definition': 'sum(stoichiometric_coefficient * energy); products positive, reactants negative',
            'energy_hartree': energy, 'energy_ev': energy * HARTREE_EV, 'conditions': reference, 'terms': rows,
            'interpretation': '음수는 선택한 에너지 정의에서 생성물 쪽 에너지가 낮음을 뜻합니다.',
            'limitations': ['전자 에너지이며 진동 ZPE/열역학 보정을 자동 포함하지 않음',
                            '촉매 활성·선택성·생산 효율 또는 반응 장벽의 증거가 아님',
                            'SCF 수렴과 기저/격자/스핀 바닥상태 검증은 별개',
                            '주기계는 동일 셀·k점·smearing 기준만 허용; 다른 참조 셀은 현재 미지원']}


def adsorption_energy(combined, substrate, adsorbate, factor=1., energy_field='total_energy_hartree'):
    if isinstance(factor, bool) or not math.isfinite(factor) or factor <= 0:
        raise ValueError('흡착종 기준 계수는 양의 유한한 수여야 합니다.')
    result = reaction_energy([(1., combined), (-1., substrate), (-factor, adsorbate)], energy_field)
    result.update(analysis='adsorption_energy',
                  definition=f'E_ads = E(combined) - E(substrate) - {factor:g} * E(adsorbate reference)',
                  interpretation='음수는 선택한 기준에서 전자적으로 유리한 흡착입니다. 촉매 활성의 증거는 아닙니다.')
    return result


def compare_runs(paths):
    if len(paths) < 2:
        raise ValueError('비교할 기록을 두 개 이상 선택하세요.')
    rows = []
    for path in paths:
        result, source = read_result(path)
        rows.append({'formula': result['formula'], 'engine': result.get('engine'), 'status': result['status'],
                     'energy_hartree': result.get('total_energy_hartree'), 'config': result['config'], 'source': source})
    return {'analysis': 'run_comparison', 'runs': rows,
            'note': '조건별 원시 결과 비교입니다. 다른 구조/방법의 에너지 차이를 반응·흡착에너지로 해석하지 마세요.'}


def save_analysis(result, root='outputs/analysis'):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    directory = Path(mkdtemp(prefix=datetime.now().strftime('%Y%m%d_%H%M%S_'), dir=root))
    path = directory / 'analysis.json'
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    return path


def main():
    parser = argparse.ArgumentParser(description='조건·원소·전하 수지를 검사하는 DFT 에너지 분석')
    parser.add_argument('--combined', type=Path)
    parser.add_argument('--substrate', type=Path)
    parser.add_argument('--adsorbate', type=Path)
    parser.add_argument('--factor', type=float, default=1.)
    parser.add_argument('--terms', type=Path, help='[{"coefficient": 1, "path": "run"}, ...] JSON')
    parser.add_argument('--energy-field', choices=ENERGY_FIELDS, default=ENERGY_FIELDS[0])
    parser.add_argument('--output', type=Path, default=Path('outputs/analysis'))
    args = parser.parse_args()
    try:
        if args.terms:
            if args.combined or args.substrate or args.adsorbate:
                parser.error('--terms와 흡착에너지 옵션을 함께 쓰지 마세요.')
            terms = json.loads(args.terms.read_text(encoding='utf-8-sig'))
            result = reaction_energy([(item['coefficient'], item['path']) for item in terms], args.energy_field)
        else:
            if not all((args.combined, args.substrate, args.adsorbate)):
                parser.error('--combined / --substrate / --adsorbate가 필요합니다.')
            result = adsorption_energy(args.combined, args.substrate, args.adsorbate, args.factor, args.energy_field)
        path = save_analysis(result, args.output)
        print(f"{result['definition']}\n{result['energy_ev']:.9f} eV\n{result['interpretation']}\n저장: {path}")
        return 0
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print('분석 오류:', exc)
        return 1


if __name__ == '__main__':
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, 'reconfigure'):
            stream.reconfigure(encoding='utf-8', errors='replace')
    raise SystemExit(main())
