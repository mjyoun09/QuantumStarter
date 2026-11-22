"""Frozen ideal-gas mixtures from QuantumStarter's exported one-bar NASA7 species."""
import argparse
import json
import math
from pathlib import Path

from matplotlib.figure import Figure
from scipy.constants import R

from engines.gas_thermo import T0, T1, TM, nasa_values, write_csv
from workflow import new_run_dir

ROOT = Path(__file__).resolve().parent
DEFAULT_DATA = ROOT / 'docs/validation/nasa7.yaml'
DEFAULT_MIX = 'N2=0.79,O2=0.21'
# Conventional atomic weights (g/mol), as used by Cantera for natural elements.
WEIGHTS = {'C': 12.011, 'N': 14.007, 'O': 15.999}
FORMULAS = {'N2': {'N': 2}, 'O2': {'O': 2}, 'CO2': {'C': 1, 'O': 2}}
TEMPERATURES = (T0, *range(300, 2001, 50))


def parse_fractions(text):
    fractions = {}
    for pair in text.split(','):
        name, sep, value = pair.partition('=')
        name = name.strip().upper()
        if not sep or name not in FORMULAS or name in fractions:
            raise ValueError('각 기체를 한 번만 지정하세요. 예: N2=0.79,O2=0.21')
        number = float(value)
        if not math.isfinite(number) or not 0 <= number <= 1:
            raise ValueError('몰분율은 0~1 사이의 유한한 수입니다.')
        fractions[name] = number
    total = math.fsum(fractions.values())
    if abs(total - 1) > 1e-8:
        raise ValueError(f'몰분율 합은 1이어야 합니다. 입력 합 {total:g}')
    return {name: x / total for name, x in fractions.items() if x > 0}


def load_species(paths):
    entries = {}
    for path in paths:
        for item in json.loads(Path(path).read_text(encoding='utf-8-sig'))['species']:
            name, thermo = item['name'], item['thermo']
            if (name not in FORMULAS or name in entries or item['composition'] != FORMULAS[name]
                    or thermo['model'] != 'NASA7' or thermo['reference-pressure'] != '1 bar'
                    or thermo['temperature-ranges'] != [T0, TM, T1]):
                raise ValueError(f'{path}: 기체 이름·중복·조성 또는 NASA7 기준 범위를 확인하세요 ({name}).')
            blocks = thermo['data']
            if len(blocks) != 2 or any(len(b) != 7 or not all(math.isfinite(v) for v in b) for b in blocks):
                raise ValueError(f'{name}: 유한한 NASA7 계수 두 세트가 필요합니다.')
            entries[name] = item
    return entries


def mixture_report(composition, files=(DEFAULT_DATA,), temperatures=TEMPERATURES):
    fractions = parse_fractions(composition)
    species = load_species(files)
    if not fractions.keys() <= species.keys():
        raise ValueError('입력 조성에 해당하는 기체의 계수가 누락됐습니다.')
    molar_mass = sum(x * sum(WEIGHTS[e]*n for e, n in FORMULAS[k].items())
                     for k, x in fractions.items()) / 1000
    rows = []
    for temperature in temperatures:
        t = float(temperature)
        if not T0 <= t <= T1:
            raise ValueError('NASA7 온도 범위는 298.15~2000 K입니다.')
        cp = sum(x * nasa_values(species[k]['thermo']['data'][int(t > TM)], t)[0]
                 for k, x in fractions.items())
        if not math.isfinite(cp) or cp <= R:
            raise ValueError(f'{t:g} K에서 정적 열용량이 양수가 아닙니다. 계수를 확인하세요.')
        cv = cp - R
        gamma = cp / cv
        speed = math.sqrt(gamma * R * t / molar_mass)
        constant = math.sqrt(1.4 * R * t / molar_mass)
        rows.append(dict(temperature_K=t, cp_J_mol_K=cp, cv_J_mol_K=cv,
                         cp_J_kg_K=cp/molar_mass, cv_J_kg_K=cv/molar_mass, gamma=gamma,
                         sound_speed_m_s=speed, sound_speed_gamma14_m_s=constant,
                         gamma14_error_percent=100*(constant/speed-1)))
    if not rows:
        raise ValueError('온도를 하나 이상 지정하세요.')
    return dict(mole_fractions=fractions, molar_mass_kg_mol=molar_mass,
                model='ideal gas, fixed mole fractions, thermally equilibrated internal modes',
                reference_pressure_Pa=100000, species=[species[k] for k in fractions], rows=rows)


def mixture_figure(report):
    fig = Figure(figsize=(9.5, 3.8), layout='constrained')
    rows = report['rows']
    ts = [r['temperature_K'] for r in rows]
    cp, gamma, speed = fig.subplots(1, 3)
    for ax, key, label in [(cp, 'cp_J_kg_K', 'Cp (J/kg/K)'), (gamma, 'gamma', 'Cp / Cv'),
                           (speed, 'sound_speed_m_s', 'Sound speed (m/s)')]:
        ax.plot(ts, [r[key] for r in rows], label='NASA7')
        ax.set(xlabel='Temperature (K)', ylabel=label)
    speed.plot(ts, [r['sound_speed_gamma14_m_s'] for r in rows], '--', label='Fixed γ = 1.4')
    speed.legend()
    fig.suptitle(' + '.join(f'{100*x:g}% {k}' for k, x in report['mole_fractions'].items()))
    return fig


def run_mixture(composition, output_root=None, files=(DEFAULT_DATA,)):
    report = mixture_report(composition, files)
    directory = new_run_dir(output_root or ROOT/'outputs/gas_mixture')
    write_csv(directory/'mixture.csv', report['rows'])
    mixture_figure(report).savefig(directory/'mixture.png', dpi=150)
    (directory/'mixture.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    return report, directory


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mix', default=DEFAULT_MIX, help='몰분율, 예: N2=0.79,O2=0.21')
    parser.add_argument('--coefficients', nargs='+', type=Path, help='QuantumStarter가 출력한 NASA7 파일')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    try:
        report, directory = run_mixture(args.mix, args.output, args.coefficients or [DEFAULT_DATA])
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.exit(1, f'혼합기체 계산 오류: {exc}\n')
    last = report['rows'][-1]
    print(f"2000 K · γ={last['gamma']:.5f} · 음속 {last['sound_speed_m_s']:.2f} m/s · "
          f"γ=1.4 상대 편차 {last['gamma14_error_percent']:+.2f}%\n{directory}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
