"""Frozen ideal-gas mixtures from QuantumStarter's exported one-bar NASA7 species."""
import argparse
import json
import math
from pathlib import Path

from matplotlib.figure import Figure
from scipy.constants import R

from engines.gas_thermo import GASES, T0, T1, TM, nasa_values, reference, write_csv
from workflow import VERSION, new_run_dir

ROOT = Path(__file__).resolve().parent
DEFAULT_DATA = ROOT / 'docs/validation/nasa7.yaml'
DEFAULT_MIX = 'N2=0.79,O2=0.21'
PRESSURE_PA = 100_000
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


def frozen_acoustics(cp, temperature, molar_mass):
    """Return molar Cv, the heat capacity ratio, and low-frequency ideal-gas sound speed."""
    if not math.isfinite(cp) or cp <= R:
        raise ValueError(f'{temperature:g} K에서 정적 열용량이 양수가 아닙니다. 계수를 확인하세요.')
    cv = cp - R
    gamma = cp / cv
    return cv, gamma, math.sqrt(gamma * R * temperature / molar_mass)


def comparison_extrema(rows, conditions):
    extrema = []
    for quantity, unit, source, difference, percent in (
        ('molar_cp', 'J/(mol K)', 'NIST Shomate Cp', 'cp_delta_J_mol_K', 'cp_delta_percent'),
        ('specific_cp', 'J/(kg K)', 'NIST Shomate Cp', 'cp_delta_J_kg_K', 'cp_delta_percent'),
        ('sound_speed', 'm/s', 'derived from NIST Cp', 'sound_speed_delta_m_s', 'sound_speed_delta_percent'),
    ):
        peak = max(rows, key=lambda r: abs(r[difference]))
        relative_peak = max(rows, key=lambda r: abs(r[percent]))
        extrema.append(dict(quantity=quantity, unit=unit, reference_type=source, **conditions,
                            max_abs_delta=abs(peak[difference]), signed_delta_at_peak=peak[difference],
                            delta_peak_temperature_K=peak['temperature_K'],
                            max_abs_percent=abs(relative_peak[percent]), signed_percent_at_peak=relative_peak[percent],
                            percent_peak_temperature_K=relative_peak['temperature_K']))
    return extrema


def mixture_report(composition, files=(DEFAULT_DATA,), temperatures=TEMPERATURES):
    fractions = parse_fractions(composition)
    species = load_species(files)
    if not fractions.keys() <= species.keys():
        raise ValueError('입력 조성에 해당하는 기체의 계수가 누락됐습니다.')
    molar_mass = sum(x * sum(WEIGHTS[e]*n for e, n in FORMULAS[k].items())
                     for k, x in fractions.items()) / 1000
    rows = []
    comparisons = []
    composition_columns = {f'x_{name}': fractions.get(name, 0.) for name in FORMULAS}
    for temperature in temperatures:
        t = float(temperature)
        if not T0 <= t <= T1:
            raise ValueError('NASA7 온도 범위는 298.15~2000 K입니다.')
        cp = sum(x * nasa_values(species[k]['thermo']['data'][int(t > TM)], t)[0]
                 for k, x in fractions.items())
        cv, gamma, speed = frozen_acoustics(cp, t, molar_mass)
        constant = math.sqrt(1.4 * R * t / molar_mass)
        rows.append(dict(temperature_K=t, cp_J_mol_K=cp, cv_J_mol_K=cv,
                         cp_J_kg_K=cp/molar_mass, cv_J_kg_K=cv/molar_mass, gamma=gamma,
                         sound_speed_m_s=speed, sound_speed_gamma14_m_s=constant,
                         gamma14_error_percent=100*(constant/speed-1)))
        ref_cp = sum(x * reference(name, t)[0] for name, x in fractions.items())
        _, ref_gamma, ref_speed = frozen_acoustics(ref_cp, t, molar_mass)
        cp_delta, speed_delta = cp - ref_cp, speed - ref_speed
        comparisons.append(dict(temperature_K=t, pressure_Pa=PRESSURE_PA, **composition_columns,
                                nasa7_cp_J_mol_K=cp, nist_cp_J_mol_K=ref_cp, cp_delta_J_mol_K=cp_delta,
                                nasa7_cp_J_kg_K=cp / molar_mass, nist_cp_J_kg_K=ref_cp / molar_mass,
                                cp_delta_J_kg_K=cp_delta / molar_mass, cp_delta_percent=100 * cp_delta / ref_cp,
                                nasa7_gamma=gamma, nist_derived_gamma=ref_gamma, nasa7_sound_speed_m_s=speed,
                                nist_derived_sound_speed_m_s=ref_speed, sound_speed_delta_m_s=speed_delta,
                                sound_speed_delta_percent=100 * speed_delta / ref_speed))
    if not rows:
        raise ValueError('온도를 하나 이상 지정하세요.')
    sources = {name: {'url': f'https://webbook.nist.gov/cgi/cbook.cgi?ID={GASES[name][2]}&Mask=1#Thermo-Gas',
                      'shomate_switch_K': GASES[name][3]} for name in fractions}
    sampled = [r['temperature_K'] for r in comparisons]
    sample_range = [min(sampled), max(sampled)]
    conditions = dict(pressure_Pa=PRESSURE_PA, **composition_columns, valid_min_K=T0, valid_max_K=T1,
                      sampled_min_K=sample_range[0], sampled_max_K=sample_range[1], sampled_points=len(sampled))
    comparison = dict(cp_source='NIST Chemistry WebBook SRD 69, gas-phase Shomate Cp (Chase, 1998)',
                      species_sources=sources, speed_source='derived from NIST Cp at fixed composition in the low-frequency ideal-gas limit',
                      difference='input NASA7 minus reference', percent_denominator='reference value',
                      valid_temperature_range_K=[T0, T1], sampled_temperature_range_K=sample_range,
                      reference_pressure_Pa=PRESSURE_PA, sampled_points=len(sampled),
                      extrema=comparison_extrema(comparisons, conditions), rows=comparisons)
    return dict(version=VERSION, mole_fractions=fractions, molar_mass_kg_mol=molar_mass, reference_comparison=comparison,
                model='ideal gas, fixed mole fractions, thermally equilibrated internal modes',
                reference_pressure_Pa=PRESSURE_PA, species=[species[k] for k in fractions], rows=rows)


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


def reference_figure(report):
    """Heat capacity and derived sound speed with separate panels for their signed residuals."""
    fig = Figure(figsize=(9.3, 5.8), layout='constrained')
    axes = fig.subplots(2, 2, sharex='col', gridspec_kw={'height_ratios': [2, 1]})
    rows = sorted(report['reference_comparison']['rows'], key=lambda row: row['temperature_K'])
    temperatures = [r['temperature_K'] for r in rows]
    for column, (model, benchmark, delta, ylabel, unit, ref_label) in enumerate((
        ('nasa7_cp_J_kg_K', 'nist_cp_J_kg_K', 'cp_delta_J_kg_K', 'Specific heat, Cp', 'J/(kg K)', 'NIST Shomate'),
        ('nasa7_sound_speed_m_s', 'nist_derived_sound_speed_m_s', 'sound_speed_delta_m_s',
         'Sound speed', 'm/s', 'From NIST Cp'),
    )):
        upper, lower = axes[:, column]
        upper.plot(temperatures, [r[model] for r in rows], color='#176389', lw=2, label='Input NASA7')
        upper.plot(temperatures, [r[benchmark] for r in rows], '--', color='#bd5625', label=ref_label)
        upper.set_ylabel(f'{ylabel} [{unit}]')
        upper.legend(frameon=False)
        lower.axhline(0, color='.55', linewidth=.8)
        lower.plot(temperatures, [r[delta] for r in rows], color='#176389', marker='.', ms=3)
        lower.set(xlabel='Temperature (K)', ylabel=f'NASA7 − reference\n[{unit}]')
    label = ' + '.join(f'{100*x:g}% {name}' for name, x in report['mole_fractions'].items())
    fig.suptitle(f'{label}  |  mole fractions, 1 bar')
    return fig


def run_mixture(composition, output_root=None, files=(DEFAULT_DATA,)):
    report = mixture_report(composition, files)
    directory = new_run_dir(output_root or ROOT/'outputs/gas_mixture')
    write_csv(directory/'mixture.csv', report['rows'])
    mixture_figure(report).savefig(directory/'mixture.png', dpi=150)
    comparison = report['reference_comparison']
    write_csv(directory/'reference_comparison.csv', comparison['rows'])
    write_csv(directory/'reference_extrema.csv', comparison['extrema'])
    reference_figure(report).savefig(directory/'reference_comparison.png', dpi=150)
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
          f"γ=1.4 상대 편차 {last['gamma14_error_percent']:+.2f}%")
    print('NIST 종별 Cp와 유도 음속 기준 · 최대 편차와 온도 (NASA7 − 기준)')
    labels = {'molar_cp': '몰 정압 열용량', 'specific_cp': '질량 정압 비열', 'sound_speed': '음속'}
    for metric in report['reference_comparison']['extrema']:
        print(f"{labels[metric['quantity']]}  {metric['signed_delta_at_peak']:+.4f} {metric['unit']} "
              f"({metric['delta_peak_temperature_K']:g} K), {metric['signed_percent_at_peak']:+.4f}% "
              f"({metric['percent_peak_temperature_K']:g} K)")
    print(directory)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
