"""Exercise the molecular dispatcher and compare spectroscopy/scan results with independent references."""
import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
import sys
import numpy as np
from pyscf import dft, lib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular import DFTConfig, parse_coordinates
from engines import calculate
from engines.pyscf_engine import HARTREE_TO_EV, build_molecule, configure_grids
from molecular_plotting import create_molecular_figure

GEOMETRIES = {
    'water':'O 0 0 0; H 0 -.757 .587; H 0 .757 .587',
    'formaldehyde':'C 0 0 0; O 0 0 1.21; H 0 .94 -.59; H 0 -.94 -.59',
    'ethylene':'C -.67 0 0; C .67 0 0; H -1.23 .93 0; H -1.23 -.93 0; H 1.23 .93 0; H 1.23 -.93 0'}


def check_spectrum(cfg, out):
    summary = calculate(cfg, out, lambda _: None)
    peaks = summary['uvvis']['peaks']
    mol = build_molecule(cfg)
    ref = dft.RKS(mol)
    configure_grids(ref.grids, cfg)
    ref.xc, ref.conv_tol, ref.conv_tol_grad, ref.small_rho_cutoff = cfg.xc, cfg.conv_tol, 1e-9, 0.
    ref.kernel()
    td = ref.TDA()
    # A full independent reference eigensystem avoids dependence on the initial
    # guess, root locking or residual tolerance of an iterative reference solver.
    a, _ = td.get_ab()
    no, nv = a.shape[:2]
    values, vectors = np.linalg.eigh(a.reshape(no*nv, no*nv))
    energies = values[:cfg.excited_states]
    xy = [(v.reshape(no, nv)/np.sqrt(2), 0) for v in vectors[:, :len(energies)].T]
    ref_f = td.oscillator_strength(e=energies, xy=xy)
    assert ref.converged and len(energies) == len(peaks)
    energy_delta = np.asarray([p['energy_eV'] for p in peaks])-energies*HARTREE_TO_EV
    strength_delta = np.asarray([p['oscillator_strength'] for p in peaks])-ref_f
    err_e, err_f = max(abs(energy_delta)), max(abs(strength_delta))
    assert err_e < 2e-5 and err_f < 2e-5, (err_e, err_f)
    create_molecular_figure(summary, out).savefig(out/'spectrum.png', dpi=150)
    return dict(states=len(peaks), max_energy_error_eV=err_e, max_strength_error=err_f,
                reference_energy_eV=(energies*HARTREE_TO_EV).tolist(), reference_strength=ref_f.tolist())


def check_gas_files(directory):
    """Read the exported species and integrate Cp/R independently of the fitter."""
    from numpy.polynomial import Polynomial
    from scipy.constants import R

    def evaluate(coeff, t):
        cp = Polynomial(coeff[:5])
        h = cp.integ(k=coeff[5])(t)
        s = coeff[0]*np.log(t) + Polynomial(coeff[1:5]).integ(k=coeff[6])(t)
        return R*np.array([cp(t), h/1000, s])

    species = json.loads((directory/'nasa7.yaml').read_text(encoding='utf-8'))['species']
    anchors = {'N2':0., 'O2':0., 'CO2':-393.5224}
    assert sorted(s['name'] for s in species) == sorted(anchors)
    with (directory/'gas_thermo.csv').open(encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 3*36
    results = {}
    for entry in species:
        name, thermo = entry['name'], entry['thermo']
        assert thermo['model'] == 'NASA7' and thermo['reference-pressure'] == '1 bar'
        assert thermo['temperature-ranges'] == [298.15, 1000., 2000.]
        low, high = thermo['data']
        assert len(low) == len(high) == 7
        assert np.max(abs(evaluate(low, 1000.)-evaluate(high, 1000.))) < 1e-8
        assert abs(evaluate(low, 298.15)[1]-anchors[name]) < 1e-8
        errors = []
        temperatures = []
        for row in (r for r in rows if r['species'] == name):
            t = float(row['temperature_K'])
            temperatures.append(t)
            fitted = evaluate(low if t <= 1000. else high, t)-[0., anchors[name], 0.]
            fields = ('cp_J_mol_K', 'delta_h_kJ_mol', 's_J_mol_K')
            tabulated = np.array([float(row['nasa7_'+k]) for k in fields])
            assert np.max(abs(fitted-tabulated)) < 1e-8
            errors.append(abs(fitted-[float(row[k]) for k in fields]))
        np.testing.assert_array_equal(temperatures, np.r_[298.15, np.arange(300., 2001., 50.)])
        maximum = np.max(errors, axis=0)
        assert np.all(maximum < [.05, .005, .02]), (name, maximum)
        results[name] = maximum.tolist()
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--out', type=Path, help='new spectroscopy and stability output directory')
    mode.add_argument('--check-gases', type=Path, help='read and verify an existing gas batch export')
    args = parser.parse_args()
    if args.check_gases:
        print(json.dumps(check_gas_files(args.check_gases), indent=2))
        return
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    lib.num_threads(1)
    base = DFTConfig(engine='native', basis='def2-svp', grid_level=3, threads=1, cube_points=12,
                     max_cycle=300, conv_tol=1e-11, native_density_tol=1e-8, native_residual_tol=1e-8)
    results = {}
    for name, xc in [('water','pbe'), ('formaldehyde','pbe'), ('ethylene','pbe'), ('water','pbe0')]:
        cfg = replace(base, atoms=parse_coordinates(GEOMETRIES[name]), xc=xc, calculation='uvvis')
        key = name+'_'+xc
        results[key] = check_spectrum(cfg, out/key)
        print(key, results[key]['max_energy_error_eV'], results[key]['max_strength_error'], flush=True)
    with (out/'uvvis_comparison.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['case', 'state', 'energy_eV', 'pyscf_energy_eV', 'wavelength_nm',
                         'oscillator_strength', 'pyscf_oscillator_strength', 'leading_transitions (MO:MO:weight)'])
        for key, data in results.items():
            with (out/key/'uvvis.csv').open(encoding='utf-8') as peaks:
                for row, energy, strength in zip(csv.DictReader(peaks), data['reference_energy_eV'], data['reference_strength']):
                    writer.writerow([key, row['state'], row['energy_eV'], energy, row['wavelength_nm'],
                                     row['oscillator_strength'], strength, row['leading_transitions (MO:MO:weight)']])
    for name, geometry, spec in [('h2','H 0 0 -.37; H 0 0 .37','distance 1 2 1.5 2.0 6'),
                                 ('water',GEOMETRIES['water'],'distance 1 2 0.95 2.3 4')]:
        cfg = replace(base, atoms=parse_coordinates(geometry), calculation='spin_scan', scan_spec=spec)
        path = out/(name+'_scan')
        summary = calculate(cfg, path, lambda _: None)
        scan = summary['spin_scan']
        assert summary['converged'] and len(scan['points']) == int(spec.split()[-1])
        assert scan['points'][0]['spin_curvature_hartree'] > 0 and scan['points'][-1]['delta_hartree'] < -1e-5
        assert len(scan['brackets_angstrom']) == 1 and scan['points'][-1]['selected'] == 'uks'
        a, b = scan['brackets_angstrom'][0]
        assert 0 < b-a <= .0004
        if name == 'h2':
            assert a < 1.61784 < b
            with (ROOT/'docs/validation/h2_scan.csv').open() as stream:
                old = {round(float(r['distance_angstrom']), 3):r for r in csv.DictReader(stream)}
            for point in scan['points']:
                ref = old.get(round(point['distance_angstrom'], 3))
                if ref:
                    assert abs(point['lowest_found_hartree']-float(ref['pyscf_uks_hartree'])) < 1e-8
        create_molecular_figure(summary, path).savefig(path/'spin_scan.png', dpi=150)
        results[name+'_scan'] = scan
        print(name, scan['brackets_angstrom'], flush=True)
    (out/'report.json').write_text(json.dumps(results, indent=2, allow_nan=False), encoding='utf-8')


if __name__ == '__main__':
    main()
