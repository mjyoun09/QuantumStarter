"""Validate radial_atom.py (log-grid spherical LDA/LSDA, H-Kr) against NIST tables and PySCF RKS.

Run in WSL: .venv-dft/bin/python scripts/validate_atoms.py [--label NAME] [--max-z 36] [--no-pyscf]

1. Every element Z = 1..max_z in LDA and LSD at grid_step 0.002 and 0.001 (default): total energy,
   its four terms (Ekin, Ecoul, Eenuc, Exc) and every orbital eigenvalue against data/nist_atomic_reference.json
   (parsed from the NIST pages by scripts/build_nist_atomic_reference.py). Richardson extrapolation
   X(0.001) + [X(0.001) - X(0.002)]/3 removes the leading h^2 error.
2. Observed order with a third grid (0.0005) for selected cases.
3. Our spin-polarized LDA formulas against Libxc (PySCF bundle) on converged LSD densities.
4. Closed-shell atoms against PySCF RKS (LDA_X,LDA_C_VWN) with even-tempered s(+p) bases.
Criteria (H-Ar): all SCF converged, default |dEtot| <= 1e-4 Ha, Richardson |dEtot| <= 2e-6 Ha and
max |d eps| <= 2e-6 Ha versus NIST (printed to 6 decimals), PySCF |dE| <= 1e-6 Ha. K-Kr is reported.
PySCF and Libxc are used only here, never inside radial_atom.py.
"""
import argparse
import csv
from datetime import datetime
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from radial_atom import SYMBOLS, XC_FUNCTIONALS, AtomConfig, lsda_xc, run_atom  # noqa: E402

COARSE, DEFAULT, FINE = 0.002, 0.001, 0.0005
TERMS = {'Etot': None, 'Ekin': 'kinetic', 'Ecoul': 'hartree', 'Eenuc': 'electron_nuclear',
         'Exc': 'exchange_correlation'}
ORDER_CASES = (('N', True), ('Ar', False), ('Kr', True))
LIBXC_SPIN_DENSITY_THRESHOLD = 1e-15
PYSCF_ATOMS = (('He', False), ('Be', False), ('Ne', True), ('Mg', True), ('Ar', True))
CRITERIA = {'default_total_hartree': 1e-4, 'richardson_total_hartree': 2e-6,
            'richardson_eigenvalue_hartree': 2e-6, 'pyscf_total_hartree': 1e-6}


def values(result):
    out = {'Etot': result.total_energy}
    out.update({term: result.components[key] for term, key in TERMS.items() if key})
    for orbital in result.orbitals:
        out[f"eps {orbital['label']} {orbital['spin']}"] = orbital['eigenvalue_hartree']
    return out


def reference(nist_column, key):
    if not key.startswith('eps '):
        return nist_column[key]
    _, label, spin = key.split()
    entry = nist_column['eigenvalues'][label]
    return entry if spin == 'unpolarized' else entry[spin]


def even_tempered(alpha_min, beta, alpha_max):
    count = int(np.ceil(np.log(alpha_max / alpha_min) / np.log(beta))) + 1
    return alpha_min * beta ** np.arange(count)


def pyscf_rks(symbol, beta, with_p):
    from pyscf import dft, gto
    basis = [[0, [float(a), 1.0]] for a in even_tempered(0.02, beta, 1e9)]
    if with_p:
        basis += [[1, [float(a), 1.0]] for a in even_tempered(0.02, beta, 1e6)]
    mol = gto.M(atom=f'{symbol} 0 0 0', basis={symbol: basis}, charge=0, spin=0, verbose=0)
    mf = dft.RKS(mol)
    mf.xc = XC_FUNCTIONALS['pyscf_xc_string']
    mf.grids.atom_grid = {symbol: (300, 50)}
    mf.grids.prune = None
    mf.small_rho_cutoff = 0.0
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    started = time.perf_counter()
    energy = float(mf.kernel())
    return {'beta': beta, 'nao': int(mol.nao), 'converged': bool(mf.converged), 'total_energy_hartree': energy,
            'basis': f"even-tempered s 0.02..1e9{' + p 0.02..1e6' if with_p else ''}, beta {beta}",
            'smallest_overlap_eigenvalue': float(np.linalg.eigvalsh(mol.intor('int1e_ovlp'))[0]),
            'seconds': time.perf_counter() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--label', default=datetime.now().strftime('%Y%m%d_%H%M%S'))
    parser.add_argument('--max-z', type=int, default=36)
    parser.add_argument('--no-pyscf', action='store_true')
    args = parser.parse_args()
    out = ROOT / 'outputs' / 'atom_validation' / args.label
    out.mkdir(parents=True, exist_ok=False)
    log = (out / 'progress.log').open('w', encoding='utf-8')

    def say(message):
        print(message, flush=True)
        log.write(message + '\n')
        log.flush()

    nist_data = json.loads((ROOT / 'data' / 'nist_atomic_reference.json').read_text(encoding='utf-8'))
    report = {'grids': {'coarse': COARSE, 'default': DEFAULT, 'fine': FINE}, 'criteria': CRITERIA,
              'nist_source': nist_data['source'], 'xc': XC_FUNCTIONALS,
              'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                              'platform': platform.platform()}}
    cases, results = [], {}
    for z in range(1, args.max_z + 1):
        symbol = SYMBOLS[z - 1]
        for column, polarized in (('LDA', False), ('LSD', True)):
            started = time.perf_counter()
            coarse = run_atom(AtomConfig(symbol, spin_polarized=polarized, grid_step=COARSE))
            default = run_atom(AtomConfig(symbol, spin_polarized=polarized, grid_step=DEFAULT))
            results[(symbol, polarized)] = default
            nist = nist_data['elements'][symbol][column]
            v_coarse, v_default = values(coarse), values(default)
            entries = {}
            for key, value in v_default.items():
                richardson = value + (value - v_coarse[key]) / 3.0
                ref = reference(nist, key)
                entries[key] = {'nist': ref, 'coarse': v_coarse[key], 'default': value, 'richardson': richardson,
                                'default_minus_nist': value - ref, 'richardson_minus_nist': richardson - ref}
            eigen = [e for k, e in entries.items() if k.startswith('eps ')]
            case = {'Z': z, 'element': symbol, 'column': column, 'configuration': default.config.configuration_text,
                    'converged': coarse.converged and default.converged,
                    'iterations': [len(coarse.history), len(default.history)],
                    'grid_points_default': int(default.grid.r.size),
                    'electron_count': default.electron_count, 'values': entries,
                    'default_total_minus_nist': entries['Etot']['default_minus_nist'],
                    'richardson_total_minus_nist': entries['Etot']['richardson_minus_nist'],
                    'richardson_max_abs_eigenvalue_minus_nist': max(abs(e['richardson_minus_nist']) for e in eigen),
                    'default_max_abs_eigenvalue_minus_nist': max(abs(e['default_minus_nist']) for e in eigen),
                    'richardson_max_abs_term_minus_nist': max(abs(entries[t]['richardson_minus_nist'])
                                                              for t in TERMS if t != 'Etot'),
                    'node_check': all(o['nodes'] == o['expected_nodes'] for o in default.orbitals),
                    'seconds': time.perf_counter() - started}
            cases.append(case)
            say(f"{z:2d} {symbol:2s} {column} conv={case['converged']} it={case['iterations']} "
                f"dE(default)={case['default_total_minus_nist']:+.2e} dE(Richardson)={case['richardson_total_minus_nist']:+.2e} "
                f"max|deps|(R)={case['richardson_max_abs_eigenvalue_minus_nist']:.1e} "
                f"max|dterm|(R)={case['richardson_max_abs_term_minus_nist']:.1e} {case['seconds']:.1f}s")
    report['nist_cases'] = cases
    with (out / 'nist_comparison.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['Z', 'element', 'column', 'quantity', 'nist', 'default', 'richardson',
                         'default_minus_nist', 'richardson_minus_nist'])
        for case in cases:
            for key, entry in case['values'].items():
                writer.writerow([case['Z'], case['element'], case['column'], key, entry['nist'], entry['default'],
                                 entry['richardson'], entry['default_minus_nist'], entry['richardson_minus_nist']])

    order = []
    for symbol, polarized in ORDER_CASES:
        if SYMBOLS.index(symbol) + 1 > args.max_z:
            continue
        fine = run_atom(AtomConfig(symbol, spin_polarized=polarized, grid_step=FINE))
        case = next(c for c in cases if c['element'] == symbol and c['column'] == ('LSD' if polarized else 'LDA'))
        e_c, e_d, e_f = case['values']['Etot']['coarse'], case['values']['Etot']['default'], fine.total_energy
        entry = {'element': symbol, 'column': case['column'], 'E_0.002': e_c, 'E_0.001': e_d, 'E_0.0005': e_f,
                 'observed_order': float(np.log2((e_c - e_d) / (e_d - e_f))),
                 'richardson_fine_pair': e_f + (e_f - e_d) / 3.0,
                 'richardson_fine_pair_minus_nist': e_f + (e_f - e_d) / 3.0 - case['values']['Etot']['nist'],
                 'fine_converged': fine.converged}
        order.append(entry)
        say(f"order {symbol} {case['column']}: p={entry['observed_order']:.3f} "
            f"R(0.001,0.0005)-NIST={entry['richardson_fine_pair_minus_nist']:+.2e}")
    report['observed_order'] = order

    if not args.no_pyscf:
        from pyscf.dft import libxc
        import pyscf
        report['environment'].update({'pyscf': pyscf.__version__, 'libxc': libxc.__version__})
        libxc_checks = []
        for symbol in ('N', 'Fe', 'Cr'):
            if SYMBOLS.index(symbol) + 1 > args.max_z:
                continue
            result = results[(symbol, True)]
            densities = result.spin_radial / (4.0 * np.pi * result.grid.r ** 2)
            keep = densities.sum(axis=0) > 1e-10
            up, down = densities[0][keep], densities[1][keep]
            exc, vxc = libxc.eval_xc(XC_FUNCTIONALS['pyscf_xc_string'], (up, down), spin=1, deriv=1)[:2]
            mine = lsda_xc(up, down)
            # Libxc 7 replaces a spin density below its threshold 1e-15 bohr^-3 by 1e-15 (measured 2026-11-17:
            # with n_up = 1e-6, Libxc at n_down = 1e-16 returns exactly our value at n_down = 1e-15). Such points
            # (outer tail of e.g. Cr 3d5 4s1 LSD, where the minority spin has no 3d/4s electron) are reported
            # separately; our analytic values there agree with a one-sided finite difference.
            away = np.minimum(up, down) >= LIBXC_SPIN_DENSITY_THRESHOLD
            relative = {'eps_xc': np.abs(mine['eps_xc'] - exc) / np.abs(exc),
                        'v_xc_majority': np.abs(mine['v_xc'][0] - vxc[0][:, 0]) / np.abs(vxc[0][:, 0]),
                        'v_xc_minority': np.abs(mine['v_xc'][1] - vxc[0][:, 1]) / np.abs(vxc[0][:, 1])}
            radial_weight = (up + down) * 4.0 * np.pi * result.grid.r[keep] ** 2 * result.grid.weights[keep]
            check = {'element': symbol, 'points': int(up.size),
                     'points_below_libxc_spin_threshold': int(np.count_nonzero(~away)),
                     'electrons_at_those_points': float(np.sum(radial_weight[~away])),
                     **{f'max_relative_{key}': float(np.max(value)) for key, value in relative.items()},
                     **{f'max_relative_{key}_both_spins_above_threshold': float(np.max(value[away]))
                        for key, value in relative.items()},
                     'xc_energy_change_if_libxc_eps_used_hartree': float(np.sum((exc - mine['eps_xc'])
                                                                                * radial_weight))}
            libxc_checks.append(check)
            say(f'libxc LSD {check}')
        report['libxc_lsda_checks'] = libxc_checks
        comparisons = []
        for symbol, with_p in PYSCF_ATOMS:
            if SYMBOLS.index(symbol) + 1 > args.max_z:
                continue
            runs = [pyscf_rks(symbol, beta, with_p) for beta in (1.8, 1.6)]
            reference_run = next((r for r in reversed(runs) if r['converged']), runs[-1])
            case = next(c for c in cases if c['element'] == symbol and c['column'] == 'LDA')
            ours = case['values']['Etot']['richardson']
            entry = {'element': symbol, 'runs': runs, 'reference_beta': reference_run['beta'],
                     'reference_converged': reference_run['converged'],
                     'pyscf_total_hartree': reference_run['total_energy_hartree'],
                     'basis_change_hartree': runs[1]['total_energy_hartree'] - runs[0]['total_energy_hartree'],
                     'quantumstarter_richardson_hartree': ours,
                     'quantumstarter_default_hartree': case['values']['Etot']['default'],
                     'richardson_minus_pyscf': ours - reference_run['total_energy_hartree'],
                     'default_minus_pyscf': case['values']['Etot']['default'] - reference_run['total_energy_hartree'],
                     'nist': case['values']['Etot']['nist']}
            comparisons.append(entry)
            say(f"pyscf {symbol}: E={entry['pyscf_total_hartree']:.9f} (beta {entry['reference_beta']}, "
                f"conv {entry['reference_converged']}, basis change {entry['basis_change_hartree']:+.1e}) "
                f"ours R-pyscf={entry['richardson_minus_pyscf']:+.2e} default-pyscf={entry['default_minus_pyscf']:+.2e}")
        report['pyscf_closed_shell'] = comparisons

    light = [c for c in cases if c['Z'] <= 18]
    heavy = [c for c in cases if c['Z'] > 18]

    def worst(group, key):
        return max((abs(c[key]) for c in group), default=None)

    summary = {
        'cases': len(cases), 'all_converged': all(c['converged'] for c in cases),
        'all_node_counts_correct': all(c['node_check'] for c in cases),
        'H_Ar': {key: worst(light, key) for key in ('default_total_minus_nist', 'richardson_total_minus_nist',
                                                    'richardson_max_abs_eigenvalue_minus_nist',
                                                    'default_max_abs_eigenvalue_minus_nist',
                                                    'richardson_max_abs_term_minus_nist')},
        'K_Kr': {key: worst(heavy, key) for key in ('default_total_minus_nist', 'richardson_total_minus_nist',
                                                    'richardson_max_abs_eigenvalue_minus_nist',
                                                    'default_max_abs_eigenvalue_minus_nist',
                                                    'richardson_max_abs_term_minus_nist')},
        'K_Kr_max_relative_default_total': max((abs(c['default_total_minus_nist'] / c['values']['Etot']['nist'])
                                                for c in heavy), default=None)}
    checks = {'H_Ar_all_converged': all(c['converged'] for c in light),
              'H_Ar_default_total': summary['H_Ar']['default_total_minus_nist'] <= CRITERIA['default_total_hartree'],
              'H_Ar_richardson_total': summary['H_Ar']['richardson_total_minus_nist']
              <= CRITERIA['richardson_total_hartree'],
              'H_Ar_richardson_eigenvalues': summary['H_Ar']['richardson_max_abs_eigenvalue_minus_nist']
              <= CRITERIA['richardson_eigenvalue_hartree']}
    if not args.no_pyscf:
        checks['pyscf_closed_shell'] = all(abs(c['richardson_minus_pyscf']) <= CRITERIA['pyscf_total_hartree']
                                           for c in report['pyscf_closed_shell'])
    summary['checks'] = checks
    summary['passed'] = all(checks.values())
    report['summary'] = summary
    (out / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=1), encoding='utf-8')
    say(json.dumps(summary, indent=1))
    say(f'report: {out / "report.json"}')
    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
