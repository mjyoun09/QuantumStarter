"""Reproducible atom PBE validation; external libraries are oracle-only.

WSL, one thread: python scripts/validate_atomic_pbe.py {lda,pbe,reference} --label NAME
Each mode writes a separate report and flushes progress after every atom.
"""
import argparse
import inspect
import importlib.util
import json
from pathlib import Path
import platform
import sys
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from radial_atom import AtomConfig, SYMBOLS, run_atom, save_atom_result
from scripts.validate_atoms import even_tempered, values

CLOSED = ('He', 'Be', 'Ne', 'Mg', 'Ar', 'Zn', 'Kr')
ORDER = (('He', False), ('Ar', False), ('Kr', False), ('H', True), ('N', True), ('Fe', True))


def pyscf_reference(symbol, beta, radial_points=600, overlap_threshold=1e-10, fractional=False):
    from pyscf import dft, gto
    from pyscf.scf import addons
    if fractional and symbol != 'B':
        raise ValueError('The equal-m occupation check is defined for B 2p1.')
    max_l = 0 if symbol in ('He', 'Be') else 2 if symbol in ('Zn', 'Kr') else 1
    z = gto.charge(symbol)
    limits = [1e5 * z * z] + [100 * z * z] * max_l
    basis = [[ell, [float(a), 1.0]] for ell in range(max_l + 1)
             for a in even_tempered(0.02, beta, limits[ell])]
    mol = gto.M(atom=f'{symbol} 0 0 0', basis={symbol: basis}, spin=int(fractional), verbose=0)
    mf = addons.frac_occ_(dft.UKS(mol)) if fractional else dft.RKS(mol)
    mf.xc = 'PBE'
    mf.grids.atom_grid = {symbol: (radial_points, 86)}
    mf.grids.prune = None
    mf.small_rho_cutoff = 0
    mf.conv_tol = 1e-10
    # Tight primitives create very high virtual energies; record both the raw
    # gradient and its gap-weighted energy estimate below. Acceptance still
    # requires independent basis convergence and 1e-6-Ha agreement.
    mf.conv_tol_grad = 1e-4
    mf.max_cycle = 150
    def orthogonalizer(s, verbose=None):
        w, u = np.linalg.eigh(s)
        return u[:, w > overlap_threshold] / np.sqrt(w[w > overlap_threshold])
    # Use the SAME space in eig and DIIS. Do not change global PySCF settings.
    mf.check_linear_dependency = orthogonalizer
    diagnostics = {}
    if fractional:
        directions = np.vstack((np.eye(3), np.array([[1, 1, 1], [-1, 1, 0]], dtype=float)))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        radii = np.geomspace(.002, 8, 25)
        probes = (radii[:, None, None] * directions).reshape(-1, 3)
        ao = mol.eval_gto('GTOval_sph', probes)
        diagnostics['max_angular_density_spread'] = 0.
        def measure(dm):
            rho = np.einsum('gi,sij,gj->sg', ao, dm, ao, optimize=True)
            spread = float(np.ptp(rho.reshape(2, len(radii), -1), axis=2).max())
            diagnostics['max_angular_density_spread'] = max(spread, diagnostics['max_angular_density_spread'])
        mf.callback = lambda state: measure(state['dm'])
    started = time.perf_counter()
    energy = float(mf.kernel())
    gradient = mf.get_grad(mf.mo_coeff, mf.mo_occ)
    if fractional:
        measure(mf.make_rdm1())
        diagnostics['occupied_per_spin'] = [o[o > 0].tolist() for o in mf.mo_occ]
        majority = sorted(diagnostics['occupied_per_spin'][0])
        diagnostics['equal_p_shell'] = bool(
            len(majority) == 5 and np.allclose(majority, [1/3, 1/3, 1/3, 1, 1])
            and np.array_equal(mf.mo_occ[1][mf.mo_occ[1] > 0], [1, 1]))
    else:
        occupied = mf.mo_occ > 0
        gaps = mf.mo_energy[~occupied, None] - mf.mo_energy[occupied]
        diagnostics['gap_weighted_gradient_energy_estimate'] = float(0.5 * np.sum(gradient.reshape(gaps.shape) ** 2 / gaps))
    return {'beta': beta, 'energy': energy, 'converged': bool(mf.converged), 'nao': mol.nao,
            'min_overlap_eigenvalue': float(np.linalg.eigvalsh(mol.intor('int1e_ovlp'))[0]),
            'grid': [radial_points, 86], 'prune': None, 'small_rho_cutoff': 0,
            'basis': {'l_max': max_l, 'alpha_min': 0.02, 'alpha_max_by_l': limits},
            'overlap_threshold': overlap_threshold, 'retained_orbitals': mf.mo_coeff.shape[-1],
            'conv_tol': mf.conv_tol, 'conv_tol_grad': mf.conv_tol_grad,
            'gradient_norm': float(np.linalg.norm(gradient)),
            **diagnostics,
            'seconds': time.perf_counter() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['lda', 'lda-source', 'pbe', 'reference', 'compare', 'fractional'])
    parser.add_argument('--label', required=True)
    parser.add_argument('--resume', action='store_true', help='resume a report after its last completed atom')
    parser.add_argument('--heavy', action='store_true', help='also run optional Zn/Kr Gaussian references')
    parser.add_argument('--radial-report', type=Path)
    parser.add_argument('--reference-report', type=Path)
    parser.add_argument('--baseline', type=Path, default=ROOT / 'outputs/atom_validation/20261117_h_kr_v3/report.json')
    args = parser.parse_args()
    if args.resume and args.mode != 'pbe':
        parser.error('--resume applies to pbe only')
    out = ROOT / 'outputs/atomic_pbe_validation' / args.label / args.mode
    out.mkdir(parents=True, exist_ok=args.resume)
    log = (out / 'progress.log').open('a' if args.resume else 'w', encoding='utf-8')
    report = {'mode': args.mode, 'environment': {'python': platform.python_version(), 'numpy': np.__version__},
              'cases': [], 'complete': False}
    if args.resume:
        report = json.loads((out / 'report.json').read_text(encoding='utf-8'))
        report['complete'] = False

    def save(message):
        print(message, flush=True)
        log.write(message + '\n')
        log.flush()
        path = out / 'report.json'
        tmp = path.with_suffix('.tmp')
        tmp.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
        tmp.replace(path)

    if args.mode == 'lda-source':
        # A historical run can use another BLAS/thread/rounding environment.
        # Compare the actual pre-change source and current engine in one process.
        source = out / 'radial_atom_26381a7.py'
        source.write_bytes(subprocess.check_output(['git', 'show', '26381a7:radial_atom.py'], cwd=ROOT))
        spec = importlib.util.spec_from_file_location('radial_atom_before_pbe', source)
        original = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = original
        spec.loader.exec_module(original)
        for symbol in SYMBOLS:
            for spin in (False, True):
                before = original.run_atom(original.AtomConfig(symbol, spin_polarized=spin))
                after = run_atom(AtomConfig(symbol, spin_polarized=spin))
                same = (values(before) == values(after) and before.history == after.history
                        and np.array_equal(before.spin_radial, after.spin_radial)
                        and all(np.array_equal(before.potentials[k], after.potentials[k]) for k in before.potentials)
                        and all(np.array_equal(before.wavefunctions[k], after.wavefunctions[k]) for k in before.wavefunctions))
                case = {'element': symbol, 'polarized': spin, 'bitwise_equal': same,
                        'converged': before.converged and after.converged,
                        'max_energy_difference': max(abs(values(before)[k] - values(after)[k]) for k in values(before))}
                report['cases'].append(case)
                save(json.dumps(case))
        report['passed'] = all(c['bitwise_equal'] and c['converged'] for c in report['cases'])
    elif args.mode == 'lda':
        baseline = json.loads(args.baseline.read_text(encoding='utf-8'))
        report['baseline'] = str(args.baseline)
        for old in baseline['nist_cases']:
            started = time.perf_counter()
            result = run_atom(AtomConfig(old['element'], spin_polarized=old['column'] == 'LSD'))
            current = values(result)
            mismatches = {key: [value, old['values'][key]['default']] for key, value in current.items()
                          if value != old['values'][key]['default']}
            case = {'element': old['element'], 'column': old['column'], 'converged': result.converged,
                    'bitwise_equal_recorded_floats': not mismatches, 'mismatches': mismatches,
                    'total_minus_nist': result.total_energy - old['values']['Etot']['nist'],
                    'electron_count': result.electron_count, 'seconds': time.perf_counter() - started}
            report['cases'].append(case)
            save(f"{case['element']} {case['column']} converged={result.converged} bitwise={not mismatches}")
        report['passed'] = len(report['cases']) == 72 and all(c['converged'] and c['bitwise_equal_recorded_floats']
                                                           for c in report['cases'])
    elif args.mode == 'pbe':
        for symbol in SYMBOLS:
            for spin in (False, True):
                if any(c['element'] == symbol and c['polarized'] == spin for c in report['cases']):
                    continue
                started = time.perf_counter()
                result = run_atom(AtomConfig(symbol, spin_polarized=spin, functional='pbe'))
                save_atom_result(result, out / f'{symbol}_{int(spin)}', time.perf_counter() - started)
                case = {'element': symbol, 'polarized': spin, 'energy': result.total_energy,
                        'converged': result.converged, 'iterations': len(result.history),
                        'electron_error': result.electron_count - result.config.nuclear_charge,
                        'energy_expression_difference': result.components['total_direct_minus_double_counting'],
                        'nodes_correct': all(o['nodes'] == o['expected_nodes'] for o in result.orbitals)}
                if (symbol, spin) in ORDER or (symbol in CLOSED and not spin):
                    coarse = run_atom(AtomConfig(symbol, functional='pbe', spin_polarized=spin, grid_step=0.002))
                    case['coarse_energy'] = coarse.total_energy
                    case['coarse_converged'] = coarse.converged
                    case['richardson'] = (4 * result.total_energy - coarse.total_energy) / 3
                    if (symbol, spin) in ORDER:
                        fine = run_atom(AtomConfig(symbol, functional='pbe', spin_polarized=spin, grid_step=0.0005))
                        case['fine_energy'], case['fine_converged'] = fine.total_energy, fine.converged
                        ratio = (coarse.total_energy - result.total_energy) / (result.total_energy - fine.total_energy)
                        case['order'] = float(np.log2(ratio)) if ratio > 0 else None
                report['cases'].append(case)
                save(json.dumps(case))
        report['passed'] = all(c['converged'] and c.get('coarse_converged', True) and c.get('fine_converged', True)
                               and c['nodes_correct'] and abs(c['electron_error']) < 1e-9
                               and abs(c['energy_expression_difference']) < 1e-9
                               and ('order' not in c or c['order'] is not None and abs(c['order'] - 2) < 0.15)
                               for c in report['cases'])
    elif args.mode == 'fractional':
        coarse, fine = [run_atom(AtomConfig('B', functional='pbe', spin_polarized=True, grid_step=h))
                        for h in (.002, .001)]
        runs = [pyscf_reference('B', ratio, fractional=True) for ratio in (1.4, 1.3)]
        richardson = (4 * fine.total_energy - coarse.total_energy) / 3
        delta = richardson - runs[-1]['energy']
        change = runs[-1]['energy'] - runs[0]['energy']
        report['cases'].append(dict(element='B', radial_richardson=richardson, runs=runs,
                                    difference=delta, basis_change=change))
        report['passed'] = (coarse.converged and fine.converged and abs(delta) < 1e-6 and abs(change) < 1e-6
                            and all(r['converged'] and r['equal_p_shell'] and r['max_angular_density_spread'] < 1e-7
                                    for r in runs))
    elif args.mode == 'compare':
        radial_path = args.radial_report or out.parent / 'pbe/report.json'
        reference_path = args.reference_report or out.parent / 'reference/report.json'
        radial = json.loads(radial_path.read_text(encoding='utf-8'))
        reference_data = json.loads(reference_path.read_text(encoding='utf-8'))
        report.update(radial_report=str(radial_path), reference_report=str(reference_path),
                      absolute_energy_tolerance_hartree=1e-6, basis_change_tolerance_hartree=1e-6)
        for ref in reference_data['cases']:
            case = next(c for c in radial['cases'] if c['element'] == ref['element'] and not c['polarized'])
            delta = case['richardson'] - ref['runs'][-1]['energy']
            passed = (case['converged'] and case['coarse_converged'] and all(r['converged'] for r in ref['runs'])
                      and abs(delta) <= 1e-6 and abs(ref['basis_change']) <= 1e-6
                      and all(r['gap_weighted_gradient_energy_estimate'] < 1e-10 for r in ref['runs']))
            report['cases'].append({'element': ref['element'], 'radial_default': case['energy'],
                                    'radial_richardson': case['richardson'], 'pyscf': ref['runs'][-1]['energy'],
                                    'difference': delta, 'basis_change': ref['basis_change'], 'passed': bool(passed)})
        report['passed'] = (radial['complete'] and radial['passed'] and reference_data['complete'] and reference_data['passed']
                            and set(CLOSED[:5]).issubset(c['element'] for c in report['cases'])
                            and all(c['passed'] for c in report['cases']))
    else:
        import pyscf
        from pyscf.dft import libxc
        from pyscf.scf import addons
        report['environment'].update(pyscf=pyscf.__version__, libxc=libxc.__version__)
        # Preserve the actual occupation algorithm, not an assumption about its name.
        (out / 'pyscf_frac_occ_source.txt').write_text(inspect.getsource(addons.frac_occ_), encoding='utf-8')
        report['open_shell_reference'] = 'Use fractional mode for the B 2p1 UKS occupation and angular-density checks.'
        for symbol in (CLOSED if args.heavy else CLOSED[:5]):
            runs = [pyscf_reference(symbol, beta) for beta in (1.4, 1.3)]
            case = {'element': symbol, 'runs': runs, 'basis_change': runs[1]['energy'] - runs[0]['energy']}
            report['cases'].append(case)
            save(json.dumps(case))
        report['passed'] = all(all(r['converged'] for r in c['runs']) and abs(c['basis_change']) < 1e-6
                               for c in report['cases'])
    report['complete'] = True
    save(f"passed={report['passed']}")
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
