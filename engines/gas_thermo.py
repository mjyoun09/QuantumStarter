"""One-bar RRHO gas properties and two-range NASA7 fits (298.15–2000 K)."""
import argparse
import csv
from collections import Counter
from dataclasses import replace
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.constants import R
from scipy.linalg import null_space

from molecular import DFTConfig, parse_coordinates, save_dft_config
from scf_control import check_cancelled
from .frequencies import BOHR_ANGSTROM, HARTREE_TO_KJ_MOL, frequency_molecule, ideal_gas_thermo
from .optimizer import optimize_molecule
from .pyscf_engine import publish_results

T0, TM, T1 = 298.15, 1000., 2000.
# NIST SRD 69 Shomate A..H, SI units. Closed upper endpoint uses the last range below 2000 K.
GASES = {
    'N2': ('N 0 0 -.55; N 0 0 .55', 0, 'C7727379', 500.,
           (28.98641, 1.853978, -9.647459, 16.63537, .000117, -8.671914, 226.4168, 0.),
           (19.50583, 19.88705, -8.598535, 1.369784, .527601, -4.935202, 212.3900, 0.)),
    'O2': ('O 0 0 -.605; O 0 0 .605', 2, 'C7782447', 700.,
           (31.32234, -20.23531, 57.86644, -36.50624, -.007374, -8.903471, 246.7945, 0.),
           (30.03235, 8.772972, -3.988133, .788313, -.741599, -11.32468, 236.1663, 0.)),
    'CO2': ('O 0 0 -1.16; C 0 0 0; O 0 0 1.16', 0, 'C124389', 1200.,
            (24.99735, 55.18696, -33.69137, 7.948387, -.136638, -403.6075, 228.2431, -393.5224),
            (58.16639, 2.720074, -.492289, .038844, -6.447293, -425.9186, 263.6125, -393.5224))}


def reference(name, t):
    _, _, _, switch, low, high = GASES[name]
    a, b, c, d, e, f, g, h = low if t < switch else high
    x = t / 1000
    return np.array([a+b*x+c*x**2+d*x**3+e/x**2, a*x+b*x**2/2+c*x**3/3+d*x**4/4-e/x+f-h,
                     a*np.log(x)+b*x+c*x**2/2+d*x**3/3-e/(2*x**2)+g])


def nasa_values(coeff, t):
    """Cp [J/mol/K], H [kJ/mol], S [J/mol/K] for one NASA7 interval."""
    n = np.arange(5)
    terms = np.asarray(coeff[:5]) * t**n
    cp = R * sum(terms)
    h = R * (t * sum(terms/(n+1)) + coeff[5]) / 1000
    s = R * (coeff[0]*np.log(t) + sum(terms[1:]/n[1:]) + coeff[6])
    return np.array([cp, h, s])


def fit_nasa(temperatures, cp, h298, s298):
    # Fit Cp/R in scaled T. Restrict the coefficient space to equal Cp at TM.
    x = np.asarray(temperatures) / TM
    powers = x[:, None] ** np.arange(5)
    design = np.hstack((powers * (x <= 1)[:, None], powers * (x > 1)[:, None]))
    space = null_space([[1.]*5 + [-1.]*5])
    p = space @ np.linalg.lstsq(design @ space, np.asarray(cp)/R, rcond=None)[0]
    blocks = []
    for i, anchor in enumerate((T0, TM)):
        coeff = [*(p[5*i:5*i+5] / TM**np.arange(5)), 0., 0.]
        if blocks:
            _, h298, s298 = nasa_values(blocks[0], anchor)
        _, h_base, s_base = nasa_values(coeff, anchor)
        coeff[5] = (h298-h_base)*1000/R
        coeff[6] = (s298-s_base)/R
        blocks.append(coeff)
    return blocks


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def gas_properties(config, directory, progress=print):
    started = perf_counter()
    directory = Path(directory)
    composition = dict(Counter(a.element for a in config.atoms))
    name = next((k for k, v in GASES.items() if Counter(a.element for a in parse_coordinates(v[0])) == composition), None)
    if name is None or config.charge or config.spin != GASES[name][1]:
        raise ValueError('중성 N2(spin=0), O2(spin=2), CO2(spin=0)를 선택하세요.')
    if config.symmetry_number not in (0, 2):
        raise ValueError('세 기체의 회전 대칭수는 2입니다.')
    directory.mkdir(parents=True, exist_ok=True)
    if (directory/'restart_input.npz').exists():
        raise ValueError('기체 열물성은 평형 구조부터 새로 계산합니다. 새 실행을 시작하세요.')
    save_dft_config(config, directory/'settings.json')

    def report(message):
        check_cancelled(directory)
        progress(message)

    opt = optimize_molecule(replace(config, calculation='optimize'), directory/'optimize', report)
    if not opt['converged']:
        raise RuntimeError('기체 평형 구조가 수렴하지 않았습니다.')
    final = DFTConfig.from_dict(opt['final_geometry_config'])
    summary = frequency_molecule(replace(final, calculation='frequencies', thermo_temperature_k=T0), directory,
                                 report, publish=False)
    freq = summary['frequencies']
    modes = freq['frequencies_cm1']
    if (not summary['converged'] or not freq['linear'] or len(modes) != 3*len(config.atoms)-5 or min(modes) <= 0
            or freq['ideal_gas_thermochemistry']['symmetry_number'] != 2):
        raise RuntimeError('선형 평형 구조·실수 진동 모드·대칭수를 확인하세요.')
    coords = np.array([a.position for a in final.atoms])/BOHR_ANGSTROM
    masses = np.array(freq['masses_amu'])

    def rrho(t):
        gas = ideal_gas_thermo(0., modes, masses, coords, config.spin+1, 2, t, 100000)
        return HARTREE_TO_KJ_MOL*np.array([1000*gas['cp_hartree_per_k'], gas['enthalpy_hartree'],
                                           1000*gas['entropy_hartree_per_k']])

    temperatures = np.r_[T0, np.arange(300., T1+1, 50.)]
    model = np.array([rrho(t) for t in temperatures])
    raw_h0 = model[0, 1]
    model[:, 1] -= raw_h0
    nist = np.array([reference(name, t) for t in temperatures])
    nist[:, 1] -= nist[0, 1].copy()
    hf = GASES[name][4][-1]
    blocks = fit_nasa(temperatures, model[:, 0], hf, model[0, 2])
    species = dict(name=name, composition=composition, thermo={'model':'NASA7', 'reference-pressure':'1 bar',
                   'temperature-ranges':[T0, TM, T1], 'data':blocks})
    # JSON flow syntax is valid YAML. Keep the standard formation-enthalpy anchor in the coefficients.
    (directory/'nasa7.yaml').write_text(json.dumps({'species':[species]}, indent=2, allow_nan=False), encoding='utf-8')
    rows = []
    for t, ours, refs in zip(temperatures, model, nist):
        row = dict(species=name, temperature_K=t)
        fit = nasa_values(blocks[int(t > TM)], t) - [0, hf, 0]
        for j, key in enumerate(('cp_J_mol_K','delta_h_kJ_mol','s_J_mol_K')):
            row[key], row['nist_'+key], row['nasa7_'+key] = ours[j], refs[j], fit[j]
        rows.append(row)
    write_csv(directory/'gas_thermo.csv', rows)
    checks = (temperatures[1:]+temperatures[:-1])/2
    errors = [nasa_values(blocks[int(t > TM)], t) - rrho(t) - [0, hf-raw_h0, 0] for t in checks]
    summary['gas_thermo'] = dict(species=name, pressure_pa=100000, rows=rows, nasa7=species,
                                 formation_enthalpy_298_kJ_mol=hf, symmetry_number=2,
                                 nist_url='https://webbook.nist.gov/cgi/cbook.cgi?Mask=1&ID='+GASES[name][2],
                                 max_nist_errors=np.max(abs(model-nist),axis=0).tolist(),
                                 max_nasa_midpoint_errors=np.max(np.abs(errors),axis=0).tolist())
    summary['config'], summary['calculation'] = config.to_dict(), 'gas_thermo'
    summary['final_geometry_config'] = final.to_dict()
    summary['assumptions'][1] = 'one-bar ideal-gas RRHO at the optimized geometry'
    summary['elapsed_seconds'] = perf_counter() - started
    publish_results(summary, directory)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='new directory')
    out = parser.parse_args().out
    out.mkdir(parents=True, exist_ok=False)
    tables, species = [], []
    for name, (coords, spin, *_) in GASES.items():
        cfg = DFTConfig(atoms=parse_coordinates(coords), spin=spin, engine='native', calculation='gas_thermo',
                        xc='pbe', basis='def2-tzvp', opt_convergence='tight', grid_level=4, cube_points=12,
                        threads=1, conv_tol=1e-10, native_density_tol=1e-8, native_residual_tol=1e-8,
                        thermo_pressure_pa=100000)
        block = gas_properties(cfg, out/name)['gas_thermo']
        tables.extend(block['rows'])
        species.append(block['nasa7'])
    write_csv(out/'gas_thermo.csv', tables)
    (out/'nasa7.yaml').write_text(json.dumps({'species':species}, indent=2, allow_nan=False), encoding='utf-8')


if __name__ == '__main__':
    main()
