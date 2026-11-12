"""Independent check of the radial-grid He LDA SCF (atomic_dft.py) against PySCF RKS.

Run in WSL: .venv-dft/bin/python scripts/validate_he_lda.py [--label NAME] [--compare-result results.json ...]

1. QuantumStarter default run (the CLI default configuration).
2. Grid-spacing series at fixed r_max = 30 bohr (dr = 0.004 ... 0.000125): observed order, Richardson.
3. Domain series at fixed dr = 0.001 bohr (r_max = 5 ... 40 bohr).
4. Electron count, and a small iteration limit that must be reported as not converged.
5. Same functional: our closed-form LDA vs Libxc on the converged He density.
6. PySCF RKS with xc 'LDA_X,LDA_C_VWN' (Libxc IDs 1 and 7): even-tempered s basis series, radial and
   angular grid series, added p functions, and two standard contracted bases for context.
   Reference = the most complete even-tempered basis whose SCF converged.
7. Differences of QuantumStarter energies (default, finest grid, Richardson, optional saved runs) from
   the reference; target |dE| <= 1e-4 hartree for the default run.

No reference value is hardcoded; every number is computed here and written to report.json.
PySCF is used only in this comparison, never inside atomic_dft.py.
"""
import argparse
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from atomic_dft import XC_FUNCTIONALS, AtomicDFTConfig, lda_xc, run_scf  # noqa: E402

TARGET_ERROR_HARTREE = 1e-4
GRID_SERIES_POINTS = (7499, 14999, 29999, 59999, 119999, 239999)   # r_max 30 -> dr 0.004 ... 0.000125
DOMAIN_SERIES_R_MAX = (5.0, 7.5, 10.0, 12.5, 15.0, 20.0, 30.0, 40.0)
SERIES_TOLERANCES = {'energy_tol_hartree': 1e-12, 'density_tol_electrons': 1e-7}
# (alpha_min, beta, alpha_max) in bohr^-2; ordered from least to most complete.
S_BASIS_SERIES = ((0.05, 2.5, 1e8), (0.05, 2.0, 1e8), (0.02, 2.0, 1e8), (0.02, 1.8, 1e8), (0.02, 1.6, 1e9))


def even_tempered(alpha_min, beta, alpha_max):
    count = int(np.ceil(np.log(alpha_max / alpha_min) / np.log(beta))) + 1
    return alpha_min * beta ** np.arange(count)


def ours(config, label):
    started = time.perf_counter()
    result = run_scf(config)
    return {'label': label, 'dr_bohr': result.dr, 'r_max_bohr': config.r_max_bohr,
            'grid_points': config.grid_points, 'status': result.status, 'iterations': len(result.history),
            'total_energy_hartree': result.total_energy, 'eigenvalue_1s_hartree': result.eigenvalue_1s,
            'electron_count': result.electron_count,
            'final_density_residual_electrons': result.history[-1]['density_residual_electrons'],
            'seconds': time.perf_counter() - started}, result


def pyscf_rks(s_exponents, nrad=300, nang=50, p_exponents=None, standard_basis=None):
    from pyscf import dft, gto
    if standard_basis is None:
        basis = [[0, [float(a), 1.0]] for a in s_exponents]
        basis += [[1, [float(a), 1.0]] for a in (p_exponents if p_exponents is not None else ())]
        description = (f'even-tempered s: {len(s_exponents)} functions {s_exponents[0]:.3g}..{s_exponents[-1]:.3g}'
                       + (f' + {len(p_exponents)} p' if p_exponents is not None else ''))
    else:
        basis, description = standard_basis, standard_basis
    started = time.perf_counter()
    mol = gto.M(atom='He 0 0 0', basis={'He': basis}, charge=0, spin=0, verbose=0)
    mf = dft.RKS(mol)
    mf.xc = XC_FUNCTIONALS['pyscf_xc_string']
    mf.grids.atom_grid = {'He': (nrad, nang)}
    mf.grids.prune = None
    mf.small_rho_cutoff = 0.0
    mf.conv_tol = 1e-12
    # PySCF 2.14 handles near-linear dependence inside SCF by default (remove_linear_dep_ is deprecated);
    # the smallest overlap eigenvalue is recorded so that ill-conditioned bases can be recognized.
    energy = float(mf.kernel())
    return {'basis': description, 'nao': int(mol.nao), 'radial_points': nrad, 'angular_points': nang,
            'converged': bool(mf.converged), 'total_energy_hartree': energy,
            'homo_eigenvalue_hartree': float(mf.mo_energy[mf.mo_occ > 0][-1]),
            'electron_count_grid': float(mf.grids.weights @ mf._numint.get_rho(mol, mf.make_rdm1(), mf.grids)),
            'smallest_overlap_eigenvalue': float(np.linalg.eigvalsh(mol.intor('int1e_ovlp'))[0]),
            'seconds': time.perf_counter() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--label', default=datetime.now().strftime('%Y%m%d_%H%M%S'))
    parser.add_argument('--compare-result', type=Path, action='append', default=[],
                        help='saved he-lda results.json (e.g. the Windows CLI run) to compare as well')
    args = parser.parse_args()
    out = ROOT / 'outputs' / 'he_lda_validation' / args.label
    out.mkdir(parents=True, exist_ok=False)
    log = (out / 'progress.log').open('w', encoding='utf-8')

    def say(message):
        print(message, flush=True)
        log.write(message + '\n')
        log.flush()

    import pyscf
    from pyscf.dft import libxc
    report = {'system': 'neutral He, 1s^2, spin-unpolarized LDA (Slater + VWN5)', 'xc': XC_FUNCTIONALS,
              'target_error_hartree': TARGET_ERROR_HARTREE,
              'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                              'pyscf': pyscf.__version__, 'libxc': libxc.__version__, 'platform': platform.platform()},
              'pyscf_parse_xc': str(libxc.parse_xc(XC_FUNCTIONALS['pyscf_xc_string']))}

    # 1. default run
    default_config = AtomicDFTConfig()
    default, default_result = ours(default_config, 'default')
    report['default_config'] = asdict(default_config)
    report['quantumstarter_default'] = default
    say(f"default: {default['status']} {default['iterations']} it, E={default['total_energy_hartree']:.12f} "
        f"eps={default['eigenvalue_1s_hartree']:.12f} N={default['electron_count']:.12f}")

    # 2. grid series
    grid = []
    for points in GRID_SERIES_POINTS:
        entry, _ = ours(AtomicDFTConfig(r_max_bohr=30.0, grid_points=points, **SERIES_TOLERANCES), f'N={points}')
        if grid:
            entry['change_from_previous_hartree'] = entry['total_energy_hartree'] - grid[-1]['total_energy_hartree']
            entry['richardson_total_energy_hartree'] = (entry['total_energy_hartree']
                                                        + entry['change_from_previous_hartree'] / 3.0)
            entry['richardson_eigenvalue_1s_hartree'] = (entry['eigenvalue_1s_hartree']
                                                         + (entry['eigenvalue_1s_hartree']
                                                            - grid[-1]['eigenvalue_1s_hartree']) / 3.0)
        if len(grid) >= 2:
            entry['observed_order'] = float(np.log2(grid[-1]['change_from_previous_hartree']
                                                    / entry['change_from_previous_hartree']))
        grid.append(entry)
        say(f"grid dr={entry['dr_bohr']:.6g}: {entry['status']} E={entry['total_energy_hartree']:.12f} "
            f"order={entry.get('observed_order', float('nan')):.4f} "
            f"richardson={entry.get('richardson_total_energy_hartree', float('nan')):.12f}")
    report['grid_series'] = grid

    # 3. domain series
    domain = []
    for r_max in DOMAIN_SERIES_R_MAX:
        entry, _ = ours(AtomicDFTConfig(r_max_bohr=r_max, grid_points=int(round(r_max / 0.001)) - 1,
                                        **SERIES_TOLERANCES), f'r_max={r_max:g}')
        domain.append(entry)
        say(f"domain r_max={r_max:g}: {entry['status']} E={entry['total_energy_hartree']:.13f}")
    for entry in domain:
        entry['difference_from_largest_hartree'] = entry['total_energy_hartree'] - domain[-1]['total_energy_hartree']
    report['domain_series'] = domain

    # 4. non-convergence
    stopped, _ = ours(AtomicDFTConfig(max_iterations=3), 'max_iterations=3')
    report['small_iteration_limit'] = stopped
    say(f"max_iterations=3: {stopped['status']} residual={stopped['final_density_residual_electrons']:.3e}")

    # 5. same functional on the converged density
    rho = default_result.density[default_result.density > 1e-10]
    exc, vxc = libxc.eval_xc(XC_FUNCTIONALS['pyscf_xc_string'], rho, spin=0, deriv=1)[:2]
    mine = lda_xc(rho)
    report['libxc_on_he_density'] = {
        'points': int(rho.size), 'density_range': [float(rho.min()), float(rho.max())],
        'max_relative_difference_eps_xc': float(np.max(np.abs(mine['eps_xc'] - exc) / np.abs(exc))),
        'max_relative_difference_v_xc': float(np.max(np.abs(mine['v_xc'] - vxc[0]) / np.abs(vxc[0])))}
    say(f"libxc on He density: {report['libxc_on_he_density']}")

    # 6. PySCF
    basis_series = []
    for alpha_min, beta, alpha_max in S_BASIS_SERIES:
        entry = pyscf_rks(even_tempered(alpha_min, beta, alpha_max))
        entry['even_tempered'] = {'alpha_min': alpha_min, 'beta': beta, 'alpha_max': alpha_max}
        basis_series.append(entry)
        say(f"pyscf {entry['basis']}: conv={entry['converged']} E={entry['total_energy_hartree']:.12f} "
            f"homo={entry['homo_eigenvalue_hartree']:.12f} Smin={entry['smallest_overlap_eigenvalue']:.1e}")
    report['pyscf_basis_series'] = basis_series
    converged = [entry for entry in basis_series if entry['converged']]
    reference = converged[-1]
    reference_exponents = even_tempered(*(reference['even_tempered'][k] for k in ('alpha_min', 'beta', 'alpha_max')))
    grid_checks = [pyscf_rks(reference_exponents, nrad=nrad, nang=nang)
                   for nrad, nang in ((75, 50), (150, 50), (600, 50), (1000, 50), (300, 302))]
    grid_checks.append(pyscf_rks(reference_exponents, p_exponents=even_tempered(0.1, 2.5, 50.0)))
    for entry in grid_checks:
        say(f"pyscf check {entry['basis']} rad={entry['radial_points']} ang={entry['angular_points']}: "
            f"conv={entry['converged']} E={entry['total_energy_hartree']:.12f}")
    report['pyscf_grid_and_p_checks'] = grid_checks
    report['pyscf_standard_bases_context'] = [pyscf_rks(None, standard_basis=name)
                                              for name in ('cc-pvtz', 'cc-pv5z', 'aug-cc-pv5z')]
    for entry in report['pyscf_standard_bases_context']:
        say(f"pyscf context {entry['basis']}: E={entry['total_energy_hartree']:.12f}")
    reference_energy = reference['total_energy_hartree']
    previous = converged[-2]['total_energy_hartree'] if len(converged) > 1 else float('nan')
    report['pyscf_reference'] = {
        'basis': reference['basis'], 'total_energy_hartree': reference_energy,
        'homo_eigenvalue_hartree': reference['homo_eigenvalue_hartree'],
        'basis_change_from_previous_converged_hartree': reference_energy - previous,
        'max_grid_or_p_change_hartree': max(abs(e['total_energy_hartree'] - reference_energy) for e in grid_checks),
        'selection': 'most complete even-tempered s basis in the series whose SCF converged (conv_tol 1e-12)'}

    # 7. comparison
    finest = grid[-1]
    rows = [('default (dr=0.001, r_max=30)', default['total_energy_hartree'], default['eigenvalue_1s_hartree']),
            (f"finest grid (dr={finest['dr_bohr']:.6g})", finest['total_energy_hartree'],
             finest['eigenvalue_1s_hartree']),
            ('Richardson (last two grids)', finest['richardson_total_energy_hartree'],
             finest['richardson_eigenvalue_1s_hartree'])]
    for path in args.compare_result:
        saved = json.loads(Path(path).read_text(encoding='utf-8'))
        rows.append((f"saved run {saved['environment'].get('platform', '')} {saved['status']} ({path.name})",
                     saved['total_energy_hartree'], saved['eigenvalue_1s_hartree']))
    report['comparison'] = [{'quantumstarter': name, 'total_energy_hartree': energy,
                             'difference_from_pyscf_hartree': energy - reference_energy,
                             'eigenvalue_1s_hartree': eigenvalue,
                             'eigenvalue_difference_from_pyscf_homo_hartree':
                                 eigenvalue - reference['homo_eigenvalue_hartree']}
                            for name, energy, eigenvalue in rows]
    report['default_within_target'] = abs(default['total_energy_hartree'] - reference_energy) <= TARGET_ERROR_HARTREE
    for row in report['comparison']:
        say(f"compare {row['quantumstarter']}: E={row['total_energy_hartree']:.12f} "
            f"dE={row['difference_from_pyscf_hartree']:+.3e} deps={row['eigenvalue_difference_from_pyscf_homo_hartree']:+.3e}")
    say(f"PySCF reference E={reference_energy:.12f} ({reference['basis']}); "
        f"default within {TARGET_ERROR_HARTREE:g} Ha: {report['default_within_target']}")
    (out / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    say(f'report: {out / "report.json"}')
    return 0 if report['default_within_target'] and default['status'] == 'converged' else 1


if __name__ == '__main__':
    raise SystemExit(main())
