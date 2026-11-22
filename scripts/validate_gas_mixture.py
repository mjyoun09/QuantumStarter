"""Cross-check the released NASA7 export in Cantera, including both sides of the range switch."""
import argparse
from pathlib import Path
import sys

import cantera as ct
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gas_mixture import DEFAULT_DATA, DEFAULT_MIX, TEMPERATURES, mixture_report
from engines.gas_thermo import write_csv


def validate():
    exported = DEFAULT_DATA.read_text(encoding='utf-8')
    gas = ct.Solution(thermo='ideal-gas', species=ct.Species.list_from_yaml(exported, section='species'))
    temperatures = sorted((*TEMPERATURES, 999.99, 1000.01))
    comparisons = []
    for mix in ('N2=1', 'O2=1', 'CO2=1', DEFAULT_MIX, 'N2=0.70,O2=0.20,CO2=0.10'):
        report = mixture_report(mix, temperatures=temperatures)
        errors = []
        for row in report['rows']:
            gas.TPX = row['temperature_K'], 100000., report['mole_fractions']
            ours = [row[k] for k in ('cp_J_kg_K','cv_J_kg_K','gamma','sound_speed_m_s')]
            reference = [gas.cp_mass, gas.cv_mass, gas.cp_mass/gas.cv_mass, gas.sound_speed]
            np.testing.assert_allclose(ours, reference, rtol=2e-10, atol=1e-11)
            errors.append(np.abs(np.subtract(ours, reference)))
        cp, cv, gamma, speed = np.max(errors, axis=0)
        comparisons.append(dict(mix=mix, temperatures=len(temperatures), cantera=ct.__version__,
                                max_cp_error_J_kg_K=cp, max_cv_error_J_kg_K=cv,
                                max_gamma_error=gamma, max_speed_error_m_s=speed))
    return comparisons


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path)
    args = parser.parse_args()
    rows = validate()
    for row in rows:
        print(row)
    if args.csv:
        write_csv(args.csv, rows)
