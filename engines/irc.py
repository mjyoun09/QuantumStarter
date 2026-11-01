"""Intrinsic reaction coordinate from a transition state (geomeTRIC IRC + this project's forces).

geomeTRIC 1.1 follows the mass-weighted steepest-descent path with the two-step
method of C. Gonzalez and H. B. Schlegel, J. Phys. Chem. 94, 5523 (1990), in both
directions of the imaginary mode, and switches to an ordinary minimization near each
end. Its starting Hessian is a finite difference of this engine's analytic gradients.
geomeTRIC reports overall success even when one direction used up its step limit, so each
end is checked against the gradient criteria here; an end that fails is finished by an
ordinary minimization from that IRC structure (recorded beside the IRC end itself).
The path, the two end structures and their energies relative to the TS are recorded;
whether the ends are the intended reactant and product is for the user to compare.
"""
from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import numpy as np

from molecular import OPT_CONVERGENCE, Atom, save_dft_config
from .optimizer import _geometric_molecule, _recording_engine, _start, _xyz, optimize_molecule
from .pyscf_engine import BOHR_ANGSTROM, calculate_molecule

HARTREE_TO_KJ_MOL = 2625.4996394799


def superposition_distance(reference, moving, weights, allow_reflection=False):
    """sqrt(sum_i w_i |R x_i + t - y_i|^2) after the weighted least-squares rigid fit (Kabsch).

    With allow_reflection the fit may also invert the structure, so mirror images match."""
    w = np.asarray(weights, dtype=float)[:, None]
    p = np.asarray(reference, dtype=float)
    q = np.asarray(moving, dtype=float)
    p = p - (w * p).sum(0) / w.sum()
    q = q - (w * q).sum(0) / w.sum()
    u, _, vt = np.linalg.svd((w * q).T @ p)
    d = np.ones(3)
    if not allow_reflection and np.linalg.det(u @ vt) < 0:
        d[2] = -1.0
    return float(np.sqrt((w * (q @ (u @ np.diag(d) @ vt) - p) ** 2).sum()))


def rmsd(reference, moving, allow_reflection=False):
    """Unweighted best-fit RMSD in the input length unit."""
    n = len(reference)
    return superposition_distance(reference, moving, np.ones(n), allow_reflection) / np.sqrt(n)


def _atom_gradient_norms(gradient):
    norms = np.linalg.norm(np.asarray(gradient, dtype=float).reshape(-1, 3), axis=1)
    return float(np.sqrt(np.mean(norms**2))), float(norms.max())


def _path_from_xyz(file, trajectory, gradients, natoms):
    """Rebuild the path geomeTRIC last wrote (Å, 10 decimals) with full-precision energies and
    gradients from the recorded calculations; None if any frame has no matching calculation."""
    if not file.is_file():
        return None
    lines = file.read_text(encoding='utf-8').splitlines()
    steps = [row for row in trajectory if row['kind'] == 'step']
    xyzs, energies, grads = [], [], []
    for first in range(0, len(lines) - natoms - 1, natoms + 2):
        x = np.array([[float(v) for v in line.split()[1:4]] for line in lines[first + 2:first + 2 + natoms]])
        row = next((r for r in steps if np.abs(np.array(r['coordinates_angstrom']) - x).max() < 1e-7), None)
        if row is None:
            return None
        xyzs.append(x)
        energies.append(row['energy_hartree'])
        grads.append(gradients[row['step']])
    return (xyzs, energies, grads) if xyzs else None


def irc_molecule(config, directory, progress=print):
    try:
        import geometric
        import geometric.engine
        import geometric.molecule
        import geometric.optimize
        from geometric.errors import GeomOptNotConvergedError, IRCError
        from geometric.molecule import PeriodicTable
    except ImportError as exc:
        raise RuntimeError('IRC에는 geomeTRIC이 필요합니다. '
                           'wsl --exec .venv-dft/bin/python -m pip install -r requirements-optimize.txt') from exc
    from .frequencies import harmonic_analysis, isotope_masses
    directory, initial, step_config = _start(config, directory, 'IRC')
    criteria = OPT_CONVERGENCE[config.opt_convergence]
    elements = [a.element for a in initial]
    molecule = _geometric_molecule(geometric, initial)
    # run_optimizer checks the fragment count (M.molecules) for IRC before it builds any
    # coordinates; a file-read Molecule has it, a custom-engine one must build it first.
    molecule.build_topology()
    engine, trajectory, state = _recording_engine(geometric.engine.Engine, molecule, step_config, directory,
                                                  progress, 'IRC')
    progress(f'geomeTRIC {geometric.__version__} IRC · 질량가중 Gonzalez–Schlegel · 양 방향 · '
             f'끝점 기준 {config.opt_convergence} · 방향당 최대 {config.opt_max_steps}단계 · '
             f'시작 Hessian은 해석적 힘의 유한차분({6 * len(initial)}회)')
    converged, path = True, None
    try:
        path = geometric.optimize.run_optimizer(
            customengine=engine, input=str(directory / 'geometric'), maxiter=config.opt_max_steps,
            logIni=str(Path(__file__).with_name('geometric_log.ini')),
            convergence_energy=criteria['energy'], convergence_grms=criteria['grms'],
            convergence_gmax=criteria['gmax'], convergence_drms=criteria['drms'], convergence_dmax=criteria['dmax'],
            irc=True, irc_direction='both', hessian='first')
    except IRCError as exc:
        raise ValueError(f'IRC 시작 구조가 허수 진동 1개인 1차 안장점이 아닙니다(geomeTRIC: {str(exc).strip()}). '
                         'transition_state 계산으로 구조를 정밀화하고 진동수로 허수 진동 1개를 확인하세요.') from exc
    except GeomOptNotConvergedError:
        converged = False
    if not trajectory:
        raise RuntimeError('geomeTRIC이 계산을 한 번도 요청하지 않았습니다. geometric.log를 확인하세요.')

    # The TS Hessian geomeTRIC built from our gradients, re-analysed with this project's
    # masses: an independent count of imaginary modes at the starting structure.
    ts_coords = np.array([a.position for a in initial]) / BOHR_ANGSTROM
    hessian_files = sorted((directory / 'geometric.tmp').rglob('hessian.txt'))
    ts_modes = None
    if hessian_files:
        analysis = harmonic_analysis(np.loadtxt(hessian_files[0]), ts_coords, isotope_masses(elements))
        frequencies = np.asarray(analysis['frequencies_cm1'])
        ts_modes = {'frequencies_cm1': frequencies.tolist(),
                    'imaginary_frequencies_cm1': frequencies[frequencies < 0].tolist(),
                    'source': 'geomeTRIC finite-difference Hessian of this engine\'s analytic gradients, '
                              'analysed with QuantumStarter harmonic analysis (common-isotope masses)'}

    if path is not None:
        points = [np.asarray(x, dtype=float) for x in path.xyzs], list(path.qm_energies), list(path.qm_grads)
        path_source = 'geomeTRIC IRC path (accepted steps)'
    else:  # geomeTRIC raised (step limit in its last direction) and its path object is lost.
        points = _path_from_xyz(directory / 'geometric_irc.xyz', trajectory, state['gradients'], len(initial))
        path_source = ('geometric_irc.xyz that geomeTRIC wrote before stopping, matched to the recorded '
                       'calculations for full-precision energies and gradients') if points else None
    record = None
    if points is not None:
        xyzs, energies, gradients = points
        start = np.array([a.position for a in initial])
        matches = [i for i, x in enumerate(xyzs) if np.abs(x - start).max() < 1e-8]
        if len(matches) != 1:
            raise RuntimeError(f'geomeTRIC IRC 경로에서 TS 구조를 하나로 특정할 수 없습니다(일치 {len(matches)}개).')
        ts_index = matches[0]
        # geomeTRIC defines the mass-weighted path with its own masses (standard atomic weights).
        masses = np.array([PeriodicTable[e] for e in elements])
        s = [0.0]
        for previous, current in zip(xyzs, xyzs[1:]):
            s.append(s[-1] + superposition_distance(previous, current, masses) / BOHR_ANGSTROM)
        s = np.array(s) - s[ts_index]
        energies = np.asarray(energies, dtype=float)
        e_ts = energies[ts_index]
        frames = []
        for i, (x, energy, gradient) in enumerate(zip(xyzs, energies, gradients)):
            rms, largest = _atom_gradient_norms(gradient)
            frames.append({'frame': i, 'branch': 'ts' if i == ts_index else 'forward' if i < ts_index else 'backward',
                           's_bohr_amu05': float(s[i]), 'energy_hartree': float(energy),
                           'relative_to_ts_kj_mol': float((energy - e_ts) * HARTREE_TO_KJ_MOL),
                           'rms_atom_gradient_hartree_per_bohr': rms, 'max_atom_gradient_hartree_per_bohr': largest,
                           'coordinates_angstrom': x.tolist()})

        def uphill(branch):  # Largest energy rise while moving away from the TS; 0 for a descent.
            rises = np.diff(branch)
            return float(max(0.0, rises.max())) if len(rises) else 0.0

        def endpoint(index, name):
            row = frames[index]
            return {'frame': index, 'energy_hartree': row['energy_hartree'],
                    'relative_to_ts_kj_mol': row['relative_to_ts_kj_mol'],
                    'barrier_to_ts_kj_mol': -row['relative_to_ts_kj_mol'],
                    'path_length_bohr_amu05': abs(row['s_bohr_amu05']),
                    'rms_atom_gradient_hartree_per_bohr': row['rms_atom_gradient_hartree_per_bohr'],
                    'max_atom_gradient_hartree_per_bohr': row['max_atom_gradient_hartree_per_bohr'],
                    'meets_gradient_criteria': bool(row['max_atom_gradient_hartree_per_bohr'] < criteria['gmax']
                                                    and row['rms_atom_gradient_hartree_per_bohr'] < criteria['grms']),
                    'moved_from_ts': index != ts_index, 'xyz_file': f'irc_{name}_end.xyz'}

        ends = {'forward': endpoint(0, 'forward'), 'backward': endpoint(len(frames) - 1, 'backward')}
        final = {'forward': xyzs[0], 'backward': xyzs[-1]}
        for name, end in ends.items():
            if end['moved_from_ts'] and not end['meets_gradient_criteria']:
                # geomeTRIC's IRC can crawl in tiny steps near some ends (seen at the linear HNC
                # end) until its step limit. The end is then finished by an ordinary minimization
                # from that IRC structure; both the IRC end and the minimum are recorded.
                progress(f'IRC {name} 끝점이 기울기 기준 미달 · 그 구조에서 일반 최소화로 끝점을 마무리합니다.')
                start_atoms = tuple(Atom(e, tuple(p)) for e, p in zip(elements, final[name]))
                relaxed = optimize_molecule(
                    replace(config, atoms=start_atoms, unit='angstrom', calculation='optimize'),
                    directory / f'endpoint_{name}',
                    lambda message: None if message.startswith(('OPT_JSON', 'SCF_JSON')) else progress(message))
                coords = np.array([a['position'] for a in relaxed['final_geometry_config']['atoms']])
                # No gradient means the final SCF failed; JSON results cannot hold inf.
                rms, largest = (_atom_gradient_norms(relaxed['gradient']['gradient_hartree_per_bohr'])
                                if relaxed['gradient'] else (None, None))
                relative = (relaxed['total_energy_hartree'] - e_ts) * HARTREE_TO_KJ_MOL
                end['irc_end'] = {key: end[key] for key in (
                    'energy_hartree', 'relative_to_ts_kj_mol', 'rms_atom_gradient_hartree_per_bohr',
                    'max_atom_gradient_hartree_per_bohr', 'meets_gradient_criteria')}
                end.update({
                    'energy_hartree': relaxed['total_energy_hartree'], 'relative_to_ts_kj_mol': relative,
                    'barrier_to_ts_kj_mol': -relative, 'rms_atom_gradient_hartree_per_bohr': rms,
                    'max_atom_gradient_hartree_per_bohr': largest,
                    'meets_gradient_criteria': bool(relaxed['converged'] and largest is not None
                                                    and largest < criteria['gmax'] and rms < criteria['grms']),
                    'relaxation': {'directory': f'endpoint_{name}', 'converged': relaxed['converged'],
                                   'steps': relaxed['optimization']['steps'],
                                   'rmsd_from_irc_end_angstrom': rmsd(final[name], coords),
                                   'energy_change_hartree': relaxed['total_energy_hartree']
                                   - end['irc_end']['energy_hartree']}})
                final[name] = coords
        for name, end in ends.items():
            end['coordinates_angstrom'] = np.asarray(final[name]).tolist()
        a, b = final['forward'], final['backward']
        record = {'ts_frame': ts_index, 'frames': frames, 'endpoints': ends,
                  'reaction_energy_forward_to_backward_kj_mol':
                      ends['backward']['relative_to_ts_kj_mol'] - ends['forward']['relative_to_ts_kj_mol'],
                  # Frames run forward end -> TS -> backward end; walk each branch away from the TS.
                  'max_uphill_hartree': {'forward': uphill(energies[:ts_index + 1][::-1]),
                                         'backward': uphill(energies[ts_index:])},
                  'endpoint_rmsd_angstrom': rmsd(a, b),
                  'endpoint_rmsd_with_reflection_angstrom': rmsd(a, b, allow_reflection=True)}
        for name, end in ends.items():
            (directory / end['xyz_file']).write_text(_xyz(
                [Atom(e, tuple(p)) for e, p in zip(elements, final[name])],
                f'IRC {name} end{" (minimized)" if "relaxation" in end else ""} E={end["energy_hartree"]:.12f} Ha '
                f'rel_TS={end["relative_to_ts_kj_mol"]:.6f} kJ/mol'), encoding='utf-8')
        (directory / 'irc_path.xyz').write_text(''.join(
            _xyz([Atom(e, tuple(p)) for e, p in zip(elements, row['coordinates_angstrom'])],
                 f'frame {row["frame"]} {row["branch"]} s={row["s_bohr_amu05"]:.6f} bohr*amu^0.5 '
                 f'E={row["energy_hartree"]:.12f} Ha rel_TS={row["relative_to_ts_kj_mol"]:.6f} kJ/mol')
            for row in frames), encoding='utf-8')
        with (directory / 'irc.csv').open('w', encoding='utf-8') as stream:
            stream.write('frame,branch,s_bohr_amu05,energy_hartree,relative_to_ts_kj_mol,'
                         'rms_atom_gradient_hartree_per_bohr,max_atom_gradient_hartree_per_bohr\n')
            stream.writelines(
                f'{r["frame"]},{r["branch"]},{r["s_bohr_amu05"]:.8f},{r["energy_hartree"]:.12f},'
                f'{r["relative_to_ts_kj_mol"]:.6f},{r["rms_atom_gradient_hartree_per_bohr"]:.6e},'
                f'{r["max_atom_gradient_hartree_per_bohr"]:.6e}\n' for r in frames)

    # The run folder shows the TS itself, recomputed with the requested output settings.
    ts_rows = [row for row in trajectory if row['kind'] == 'step'
               and np.abs(np.array(row['coordinates_angstrom']) - np.array([a.position for a in initial])).max() < 1e-8]
    guess = None
    if ts_rows:
        with np.load(directory / ts_rows[0]['directory'] / 'electronic_structure.npz', allow_pickle=False) as data:
            guess = np.asarray(data['density_matrix'])
    progress('TS 구조를 요청한 출력 설정으로 다시 계산합니다.')
    final_config = replace(config, atoms=initial, unit='angstrom', calculation='gradient', opt_constraints='')
    summary = calculate_molecule(final_config, directory, progress, guess)
    save_dft_config(final_config, directory / 'final_settings.json')
    save_dft_config(config, directory / 'settings.json')  # The request.
    ends_ok = record is not None and all(end['meets_gradient_criteria'] and end['moved_from_ts']
                                         for end in record['endpoints'].values())
    summary['irc'] = {
        'method': 'mass-weighted steepest-descent path from the TS in both directions '
                  '(Gonzalez–Schlegel two-step IRC, geomeTRIC), minimization near each end',
        'optimizer': 'geomeTRIC', 'optimizer_version': geometric.__version__,
        'reference': 'C. Gonzalez, H. B. Schlegel, J. Phys. Chem. 94, 5523 (1990); '
                     'L.-P. Wang, C. Song, J. Chem. Phys. 144, 214108 (2016)',
        'coordinates': 'geomeTRIC internal coordinates (TRIC; DLC for a single fragment, see geometric.log)',
        'path_masses': 'geomeTRIC PeriodicTable (standard atomic weights); they define the mass-weighted path and s',
        'direction_names': 'forward/backward follow the sign of the TS imaginary mode eigenvector (arbitrary)',
        # Success is judged on the ends (gradient criteria, after minimization if needed), the
        # same way whichever direction geomeTRIC ran out of steps in.
        'geometric_converged': converged, 'converged': bool(ends_ok), 'path_source': path_source,
        'max_steps_per_direction': config.opt_max_steps, 'criteria_name': config.opt_convergence,
        'criteria': {'energy_change_hartree': criteria['energy'], 'rms_gradient_hartree_per_bohr': criteria['grms'],
                     'max_gradient_hartree_per_bohr': criteria['gmax'], 'rms_step_angstrom': criteria['drms'],
                     'max_step_angstrom': criteria['dmax']},
        'ts_hessian': ts_modes,
        'gradient_calculations': len(trajectory),
        'hessian_gradient_calculations': sum(row['kind'] == 'hessian' for row in trajectory),
        # geomeTRIC labels each direction's folder; a direction that used up its step limit
        # is handed over silently, so the per-direction count is kept for diagnosis.
        'step_calculations_per_direction': {
            name: sum(row['kind'] == 'step' and f'irc_{name}' in row['geometric_dirname'] for row in trajectory)
            for name in ('forward', 'backward')},
        'calculations': trajectory,
        **(record or {'frames': [], 'endpoints': None}),
        'path_file': 'irc_path.xyz' if record else None, 'table_file': 'irc.csv' if record else None,
        'log_file': 'geometric.log', 'completed_at_utc': datetime.now(timezone.utc).isoformat()}
    summary['calculation'] = 'irc'
    summary['config'] = config.to_dict()
    summary['final_geometry_config'] = final_config.to_dict()
    summary['assumptions'] = ['isolated finite molecule',
                              'steepest-descent reaction path from a first-order saddle point',
                              *summary['assumptions'][2:]]
    warnings = []
    if not converged and record:
        warnings.append(f'geomeTRIC IRC가 마지막 방향에서 {config.opt_max_steps}단계 한도로 멈췄습니다. 경로는 '
                        'geomeTRIC이 남긴 geometric_irc.xyz를 기록된 계산과 맞춰 복구했습니다.')
    elif not converged:
        warnings.append(f'IRC 미수렴: geomeTRIC이 방향당 {config.opt_max_steps}단계 안에 끝나지 않았고 경로를 '
                        '복구하지 못했습니다. calculations와 geometric.log로 진단하세요.')
    if record:
        for name, end in record['endpoints'].items():
            if not end['moved_from_ts']:
                warnings.append(f'IRC {name} 방향이 TS에서 움직이지 않았습니다.')
            elif not end['meets_gradient_criteria']:
                warnings.append(f'IRC {name} 끝점이 기울기 기준을 만족하지 않습니다' +
                                ('(끝점 최소화도 실패). ' if 'relaxation' in end else '. ') + '최소점이 아닙니다.')
            if 'relaxation' in end:
                relaxation = end['relaxation']
                warnings.append(
                    f'IRC {name} 방향은 geomeTRIC IRC가 끝점 기울기 기준에 도달하지 못해(보통 방향당 '
                    f'{config.opt_max_steps}단계 한도 도달) '
                    f'그 구조에서 일반 최소화로 끝점을 마무리했습니다: {relaxation["steps"]}단계, '
                    f'구조 이동 RMSD {relaxation["rmsd_from_irc_end_angstrom"]:.2e} Å, '
                    f'ΔE {relaxation["energy_change_hartree"]:+.2e} Ha. 경로 파일의 끝 점은 최소화 전 IRC 구조입니다.')
                if relaxation['rmsd_from_irc_end_angstrom'] > 0.1:
                    warnings.append(f'IRC {name} 끝점 최소화로 구조가 0.1 Å 넘게 바뀌었습니다. IRC가 끝점 근처까지 '
                                    '가지 못했을 수 있으므로 이 끝점의 연결은 별도로 확인하세요.')
        if record['endpoint_rmsd_with_reflection_angstrom'] < 1e-2 and abs(
                record['reaction_energy_forward_to_backward_kj_mol']) < 0.05:
            warnings.append('두 방향의 끝점이 같은 구조이거나 거울상입니다. 대칭적 재배열(예: NH₃ 반전)이 아니라면 '
                            '이 TS는 의도한 서로 다른 반응물·생성물을 잇지 않을 수 있습니다.')
    if ts_modes and len([f for f in ts_modes['imaginary_frequencies_cm1'] if f < -20]) != 1:
        warnings.append(f'시작 구조 Hessian의 유의한 허수 진동이 {len(ts_modes["imaginary_frequencies_cm1"])}개입니다'
                        '(1차 안장점은 정확히 1개).')
    warnings += ['IRC 끝점은 기울기 기준을 만족한 최소점 후보입니다. 끝점 구조로 진동수 계산(허수 없음)을 하고, '
                 '의도한 반응물·생성물과 구조를 비교해 연결을 확인하세요.',
                 'forward/backward 이름은 허수 모드 고유벡터의 부호에 따른 임의 구분입니다.']
    if config.grid_level < 3:
        warnings.insert(0, f'격자 수준 {config.grid_level}: 격자 응답 항이 빠진 힘의 오차가 끝점 수렴 기준과 '
                           '비슷하거나 큽니다. IRC에는 grid_level 3 이상을 권장합니다.')
    summary['warnings'] = warnings + summary['warnings']
    summary['converged'] = summary['irc']['converged'] and summary['status'] == 'converged'
    summary['status'] = 'converged' if summary['converged'] else 'not_converged'
    temporary = directory / 'results.json.tmp'
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(directory / 'results.json')
    shutil.rmtree(directory / 'geometric.tmp', ignore_errors=True)
    if record:
        ends = record['endpoints']
        progress(f'IRC {"수렴" if summary["converged"] else "미수렴"} · TS→forward 끝 '
                 f'{ends["forward"]["relative_to_ts_kj_mol"]:.3f} kJ/mol · TS→backward 끝 '
                 f'{ends["backward"]["relative_to_ts_kj_mol"]:.3f} kJ/mol · 경로 {len(record["frames"])}점')
    else:
        progress('IRC 미수렴 · 진단 결과 저장 완료')
    return summary
