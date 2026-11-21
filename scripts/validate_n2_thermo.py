"""Compare native PBE/RRHO gas-phase N2 with NIST at 1 bar, 298.15–2000 K."""
import argparse
import csv
from dataclasses import replace
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.constants import R

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from engines import calculate
from engines.frequencies import BOHR_ANGSTROM, HARTREE_TO_KJ_MOL, ideal_gas_thermo, vibrational_thermo
from molecular import DFTConfig, parse_coordinates

T_REF = 298.15
PRESSURE_PA = 100_000
NIST_URL = 'https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1#Thermo-Gas'
CODATA_S = 191.609  # J/(mol K), uncertainty 0.004


def shomate(temperature):
    """NIST SRD 69 nitrogen gas, Chase (1998), revised fit January 2009. SI units."""
    if not T_REF <= temperature <= 2000:
        raise ValueError('Temperature outside the study interval.')
    if temperature < 500:
        a, b, c, d, e, f, g = 28.98641, 1.853978, -9.647459, 16.63537, .000117, -8.671914, 226.4168
    else:
        a, b, c, d, e, f, g = 19.50583, 19.88705, -8.598535, 1.369784, .527601, -4.935202, 212.39
    t = temperature / 1000
    cp = a + b*t + c*t**2 + d*t**3 + e/t**2
    h = a*t + b*t**2/2 + c*t**3/3 + d*t**4/4 - e/t + f  # N2 formation enthalpy coefficient is zero
    s = a*math.log(t) + b*t + c*t**2/2 + d*t**3/3 - e/(2*t**2) + g
    return cp, h, s


def temperature_sweep(freq, coords_bohr):
    wavenumbers = freq['frequencies_cm1']
    masses = np.asarray(freq['masses_amu'])
    sigma = freq['ideal_gas_thermochemistry']['symmetry_number']
    if not freq['linear'] or sigma != 2 or len(wavenumbers) != 1 or wavenumbers[0] <= 0:
        raise ValueError('Expected a linear N2 minimum with symmetry number 2 and one real mode.')
    if freq['reference_max_gradient_hartree_per_bohr'] > 3e-5:
        raise ValueError('Gradient too large to report equilibrium thermochemistry.')

    def gas(temperature):
        # The constant electronic energy cancels from ΔH. Omit it before taking differences.
        return ideal_gas_thermo(0., wavenumbers, masses, coords_bohr, 1, sigma, temperature, PRESSURE_PA)

    h0 = gas(T_REF)['enthalpy_hartree']
    reference_h0 = shomate(T_REF)[1]
    joules_per_mole = HARTREE_TO_KJ_MOL * 1000
    rows = []
    max_derivative_error = 0.
    for temperature in (T_REF, *range(300, 2001, 50)):
        result = gas(temperature)
        cv_vib = vibrational_thermo(wavenumbers, temperature)['cv_vib_hartree_per_k']
        cp = 3.5 * R + cv_vib * joules_per_mole  # linear rotor: translation + rotation + pV
        dh = (result['enthalpy_hartree'] - h0) * HARTREE_TO_KJ_MOL
        s = result['entropy_hartree_per_k'] * joules_per_mole
        ref_cp, ref_h, ref_s = shomate(temperature)
        ref_dh = ref_h - reference_h0
        rows.append(dict(temperature_K=temperature, cp_J_mol_K=cp, nist_cp_J_mol_K=ref_cp, error_cp_J_mol_K=cp-ref_cp,
                         delta_h_kJ_mol=dh, nist_delta_h_kJ_mol=ref_dh, error_delta_h_kJ_mol=dh-ref_dh,
                         s_J_mol_K=s, nist_s_J_mol_K=ref_s, error_s_J_mol_K=s-ref_s))
        # Independent numerical checks of the analytical Cp and unit conversion at fixed pressure.
        step = .01
        low, high = gas(temperature - step), gas(temperature + step)
        dh_dt = (high['enthalpy_hartree'] - low['enthalpy_hartree']) * joules_per_mole / (2*step)
        t_ds_dt = temperature * (high['entropy_hartree_per_k'] - low['entropy_hartree_per_k']) * joules_per_mole / (2*step)
        max_derivative_error = max(max_derivative_error, abs(dh_dt - cp), abs(t_ds_dt - cp))
    if not np.isfinite([list(row.values()) for row in rows]).all() or max_derivative_error > 1e-6:
        raise ValueError('RRHO thermodynamic consistency check failed.')
    # Cross-check transcription, branch selection and units against the rounded NIST WebBook table.
    for temp, cp, dh, s in [(300, 29.12, .05, 191.8), (1000, 32.69, 21.46, 228.2), (2000, 35.98, 56.14, 252.1)]:
        ref_cp, ref_h, ref_s = shomate(temp)
        if abs(ref_cp-cp) > .0051 or abs(ref_h-dh) > .0051 or abs(ref_s-s) > .051:
            raise ValueError(f'NIST table check failed at {temp} K.')
    return rows, max_derivative_error


def plot_errors(rows, path):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure(figsize=(10.2, 3.5), layout='constrained')
    FigureCanvasAgg(figure)
    temperatures = [row['temperature_K'] for row in rows]
    panels = (('error_cp_J_mol_K', 'Heat capacity', 'J/(mol K)'),
              ('error_delta_h_kJ_mol', 'Enthalpy increment', 'kJ/mol'),
              ('error_s_J_mol_K', 'Entropy', 'J/(mol K)'))
    for axis, (key, title, units) in zip(figure.subplots(1, 3, sharex=True), panels):
        axis.axhline(0, color='0.45', lw=.8)
        axis.plot(temperatures, [row[key] for row in rows], color='#125c83', marker='.', ms=4)
        axis.set(title=title, xlabel='Temperature (K)', ylabel=units, xticks=[300, 1000, 1500, 2000])
        axis.grid(axis='y', alpha=.2)
    figure.suptitle('Nitrogen at 1 bar  |  PBE/RRHO minus NIST')
    figure.savefig(path, dpi=165)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='new output directory')
    parser.add_argument('--freq-step-bohr', type=float, default=.005, help='central-difference displacement')
    args = parser.parse_args()
    out = args.out
    config = DFTConfig(atoms=parse_coordinates('N 0 0 -.55; N 0 0 .55'), engine='native', xc='pbe',
                       basis='def2-tzvp', calculation='optimize', opt_convergence='tight', conv_tol=1e-10,
                       native_density_tol=1e-8, native_residual_tol=1e-8, grid_level=4, cube_points=12, threads=1,
                       thermo_temperature_k=T_REF, thermo_pressure_pa=PRESSURE_PA, freq_step_bohr=args.freq_step_bohr)
    out.mkdir(parents=True, exist_ok=False)
    opt = calculate(config, out / 'optimize')
    if not opt['converged']:
        raise RuntimeError('N2 optimization or final SCF did not converge.')
    optimized = DFTConfig.from_dict(opt['final_geometry_config'])
    vib = calculate(replace(optimized, calculation='frequencies'), out / 'frequencies')
    if not vib['converged']:
        raise RuntimeError('N2 vibrational calculation did not converge.')
    freq = vib['frequencies']
    coords = np.array([a.position for a in optimized.atoms])
    rows, check = temperature_sweep(freq, coords / BOHR_ANGSTROM)
    with (out / 'n2_thermo.csv').open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    plot_errors(rows, out / 'n2_errors.png')
    summary = dict(config=config.to_dict(), environment=vib['environment'], nist_url=NIST_URL,
                   bond_angstrom=float(np.linalg.norm(coords[1]-coords[0])), masses_amu=freq['masses_amu'],
                   frequency_cm1=freq['frequencies_cm1'][0], symmetry_number=freq['ideal_gas_thermochemistry']['symmetry_number'],
                   max_gradient_hartree_bohr=freq['reference_max_gradient_hartree_per_bohr'], temperatures=len(rows),
                   max_derivative_check_error_J_mol_K=check, codata_entropy_298_J_mol_K=CODATA_S,
                   entropy_minus_codata_J_mol_K=rows[0]['s_J_mol_K'] - CODATA_S,
                   raw_nist_enthalpy_298_kJ_mol=shomate(T_REF)[1])
    (out / 'report.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    print(f"Saved {len(rows)} temperatures to {out} (r={summary['bond_angstrom']:.6f} A, "
          f"harmonic frequency={summary['frequency_cm1']:.2f} cm^-1).")


if __name__ == '__main__':
    main()
