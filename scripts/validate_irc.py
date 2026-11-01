"""IRC connection check on a reaction with distinct ends: HCN <-> HNC isomerization.

Chain: optimize HCN and HNC independently -> P-RFO transition state from a bent guess ->
frequencies at the TS (exactly one imaginary mode) -> IRC both ways -> each IRC end is
matched to the independently optimized isomer by best-fit RMSD and energy -> frequencies
at both IRC ends (no imaginary mode). Every calculation is a normal recorded run.

Run in WSL: .venv-dft/bin/python scripts/validate_irc.py --label <name> [--basis def2-svp]
No reference values are hardcoded; every number is computed and written to report.json.
"""
import argparse
from dataclasses import replace
from datetime import datetime
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular import Atom, DFTConfig, parse_coordinates
from engines import calculate
from engines.irc import rmsd
from engines.pyscf_engine import HARTREE_TO_EV

HARTREE_TO_KJ_MOL = 2625.4996394799
# Atom order H, C, N everywhere so structures compare one-to-one.
HCN = 'H 0 0 -1.07; C 0 0 0; N 0 0 1.16'
HNC = 'H 0 0 2.17; C 0 0 0; N 0 0 1.17'
TS_GUESS = 'H 0 1.10 0.60; C 0 0 0; N 0 0 1.18'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', default=datetime.now().strftime('%Y%m%d_%H%M%S'))
    parser.add_argument('--basis', default='def2-svp')
    parser.add_argument('--xc', default='pbe')
    parser.add_argument('--engine', default='native')
    args = parser.parse_args()
    out = ROOT / 'outputs' / 'irc_validation' / args.label
    out.mkdir(parents=True, exist_ok=True)
    base = DFTConfig(atoms=parse_coordinates(HCN), xc=args.xc, basis=args.basis, engine=args.engine,
                     grid_level=3, cube_points=12, threads=1, opt_convergence='tight')
    log = (out / 'progress.log').open('w', encoding='utf-8')

    def run(name, cfg):
        start = time.time()
        result = calculate(cfg, out / name, lambda m: None if m.startswith(('OPT_JSON', 'SCF_JSON'))
                           else log.write(f'[{name}] {m}\n') or log.flush())
        print(f'{name}: {result["status"]} E={result["total_energy_hartree"]:.10f} Ha ({time.time() - start:.1f} s)',
              flush=True)
        return result

    def geometry(result):
        return [a['position'] for a in result['final_geometry_config']['atoms']]

    def atoms(coordinates):
        return tuple(Atom(e, tuple(c)) for e, c in zip('HCN', coordinates))

    report = {'system': 'HCN <-> HNC', 'level': f'{args.engine} {args.xc}/{args.basis}, grid 3',
              'created': datetime.now().isoformat()}
    hcn = run('opt_hcn', replace(base, calculation='optimize'))
    hnc = run('opt_hnc', replace(base, atoms=parse_coordinates(HNC), calculation='optimize'))
    ts = run('ts', replace(base, atoms=parse_coordinates(TS_GUESS), calculation='transition_state',
                           opt_convergence='default'))
    ts_freq = run('ts_freq', replace(base, atoms=atoms(geometry(ts)), calculation='frequencies'))['frequencies']
    irc = run('irc', replace(base, atoms=atoms(geometry(ts)), calculation='irc'))
    info = irc['irc']
    minima = {'HCN': (hcn['total_energy_hartree'], geometry(hcn)), 'HNC': (hnc['total_energy_hartree'], geometry(hnc))}
    ends = {}
    for direction, end in info['endpoints'].items():
        coords = end['coordinates_angstrom']  # Final end (minimized if geomeTRIC's IRC stopped short).
        raw = info['frames'][end['frame']]['coordinates_angstrom']  # Last IRC structure itself.
        matches = {name: {'rmsd_angstrom': rmsd(ref, coords), 'irc_end_rmsd_angstrom': rmsd(ref, raw),
                          'energy_difference_hartree': end['energy_hartree'] - energy}
                   for name, (energy, ref) in minima.items()}
        best = min(matches, key=lambda name: matches[name]['rmsd_angstrom'])
        freq = run(f'end_{direction}_freq', replace(base, atoms=atoms(coords), calculation='frequencies'))['frequencies']
        ends[direction] = {'matches': matches, 'assigned': best, 'endpoint': end,
                           'assigned_from_irc_end_alone': min(matches, key=lambda n: matches[n]['irc_end_rmsd_angstrom']),
                           'endpoint_frequencies_cm1': freq['frequencies_cm1'],
                           'endpoint_imaginary_cm1': freq['imaginary_frequencies_cm1'],
                           'endpoint_character': freq['stationary_point_character']}
    e_ts = irc['total_energy_hartree']
    report.update({
        'optimized': {name: {'energy_hartree': energy, 'coordinates_angstrom': coords}
                      for name, (energy, coords) in minima.items()},
        'ts': {'energy_hartree': e_ts, 'coordinates_angstrom': geometry(ts), 'search_steps': ts['optimization']['steps'],
               'frequencies_cm1': ts_freq['frequencies_cm1'], 'imaginary_cm1': ts_freq['imaginary_frequencies_cm1'],
               'character': ts_freq['stationary_point_character']},
        'irc': {k: v for k, v in info.items() if k not in ('frames', 'calculations')},
        'irc_frames': len(info['frames']), 'irc_gradient_calculations': info['gradient_calculations'],
        'ends': ends,
        'connects_hcn_and_hnc': sorted(e['assigned'] for e in ends.values()) == ['HCN', 'HNC'],
        'irc_ends_alone_connect_hcn_and_hnc':
            sorted(e['assigned_from_irc_end_alone'] for e in ends.values()) == ['HCN', 'HNC'],
        'irc_converged': irc['converged'],
        'barrier_from_hcn_kj_mol': (e_ts - minima['HCN'][0]) * HARTREE_TO_KJ_MOL,
        'barrier_from_hnc_kj_mol': (e_ts - minima['HNC'][0]) * HARTREE_TO_KJ_MOL,
        'isomerization_energy_hcn_to_hnc_kj_mol': (minima['HNC'][0] - minima['HCN'][0]) * HARTREE_TO_KJ_MOL,
        'isomerization_energy_hcn_to_hnc_ev': (minima['HNC'][0] - minima['HCN'][0]) * HARTREE_TO_EV})
    (out / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('irc_converged', 'connects_hcn_and_hnc', 'irc_ends_alone_connect_hcn_and_hnc',
                                            'barrier_from_hcn_kj_mol', 'barrier_from_hnc_kj_mol',
                                            'isomerization_energy_hcn_to_hnc_kj_mol')}, indent=1))
    print('step calculations per direction', info['step_calculations_per_direction'])
    for direction, end in ends.items():
        print(direction, end['assigned'], {k: {kk: f'{vv:.3e}' for kk, vv in v.items()} for k, v in end['matches'].items()},
              end['endpoint_character'], end['endpoint_imaginary_cm1'], 'relaxation:', end['endpoint'].get('relaxation'))
    print('TS imaginary', ts_freq['imaginary_frequencies_cm1'], '→', out)
    return 0 if report['connects_hcn_and_hnc'] and irc['converged'] else 1


if __name__ == '__main__':
    sys.exit(main())
