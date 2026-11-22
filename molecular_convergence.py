"""Compare molecular observables over a Cartesian product of bases and grids."""
import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
from time import perf_counter

from dft_workflow import ROOT, CalculationError, run_calculation
from molecular import DFTConfig, load_dft_config, read_xyz
from workflow import new_run_dir

OPERATIONS = ('uvvis', 'frequencies', 'spin_scan', 'gas_thermo')
DEFAULT_BASES = ('def2-svp', 'def2-tzvp')
DEFAULT_GRIDS = (2, 3)
RUN_FIELDS = ('basis', 'grid', 'status', 'wall_seconds', 'directory', 'error')
DELTA_FIELDS = ('axis', 'from_basis', 'from_grid', 'to_basis', 'to_grid',
                'observable', 'rank', 'unit', 'from_value', 'to_value', 'delta', 'status')


def measure(summary, operation):
    """Keys include observable, energy/frequency/distance rank (1-based) and unit."""
    values = {}
    if operation == 'uvvis':
        for peak in summary['uvvis']['peaks']:
            i = peak['state']
            values['excitation', i, 'eV'] = peak['energy_eV']
            values['strength', i, '1'] = peak['oscillator_strength']
    elif operation == 'frequencies':
        for i, frequency in enumerate(sorted(summary['frequencies']['frequencies_cm1']), 1):
            values['frequency', i, 'cm^-1'] = frequency
    elif operation == 'spin_scan':
        intervals = sorted(summary['spin_scan']['brackets_angstrom'])
        values['crossings', 0, 'count'] = len(intervals)
        for i, (left, right) in enumerate(intervals, 1):
            values['crossing_midpoint', i, 'angstrom'] = (left+right)/2
            values['bracket_width', i, 'angstrom'] = right-left
    else:
        for row in summary['gas_thermo']['rows']:
            values['Cp', row['temperature_K'], 'J/mol/K'] = row['cp_J_mol_K']
    return values


def difference(axis, previous, current):
    a, b = previous['values'], current['values']
    crossings = ('crossings', 0, 'count')
    unlike_branches = a.get(crossings) != b.get(crossings)
    rows = []
    for key in sorted(a.keys() | b.keys()):
        name, index, unit = key
        before, after = a.get(key), b.get(key)
        status = 'compared'
        if before is None or after is None:
            status = 'missing'
        elif unlike_branches and name in ('crossing_midpoint', 'bracket_width'):
            status = 'different_crossing_count'
        rows.append(dict(axis=axis, from_basis=previous['basis'], from_grid=previous['grid'],
                         to_basis=current['basis'], to_grid=current['grid'], observable=name, rank=index,
                         unit=unit, from_value=before, to_value=after,
                         delta=after-before if status == 'compared' else None, status=status))
    return rows


def save_tables(directory, report):
    for name, fields, rows in [('runs', RUN_FIELDS, report['runs']), ('differences', DELTA_FIELDS, report['differences'])]:
        path = directory/(name+'.csv')
        tmp = path.with_suffix('.tmp')
        with tmp.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fields, lineterminator='\n', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        tmp.replace(path)
    path = directory/'comparison.json'
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    tmp.replace(path)


def run_comparison(config, bases=DEFAULT_BASES, grids=DEFAULT_GRIDS, output_root=None, backend='auto',
                   progress=print, control=None):
    bases, grids = tuple(b.strip().lower() for b in bases), tuple(grids)
    if config.calculation not in OPERATIONS:
        raise ValueError('uvvis, frequencies, spin_scan 또는 gas_thermo 계산 설정을 선택하세요.')
    if not bases or not grids or len(bases)*len(grids) < 2:
        raise ValueError('서로 다른 계산 조건을 두 개 이상 지정하세요.')
    if len(set(bases)) != len(bases) or len(set(grids)) != len(grids):
        raise ValueError('기저와 격자 목록 안에 중복 값이 있습니다.')
    configs = {(b, g): replace(config, basis=b, grid_level=g) for b in bases for g in grids}
    # Currently each ECP element uses its common def2 core across supported def2 bases.
    models = {(c.resolved_relativity, c.resolved_nuclear_model, tuple(sorted(c.resolved_ecp)),
               c.resolved_grid_scheme) for c in configs.values()}
    if len(models) != 1:
        raise ValueError('요청한 기저 간에 상대론·핵 모형·ECP·격자 방식이 달라집니다. 동일 조건의 기저들을 선택하세요.')
    distinct = {(tuple(sorted(c.basis_by_element.items())), c.grid_level) for c in configs.values()}
    if len(distinct) != len(configs):
        raise ValueError('원소별 기저 설정을 적용하면 같은 계산 조건이 반복됩니다.')
    directory = new_run_dir(output_root or ROOT/'outputs/molecular_convergence').resolve()
    geometry = ('optimized at each condition' if config.calculation == 'gas_thermo'
                else 'identical bond path' if config.calculation == 'spin_scan' else 'fixed input geometry')
    report = dict(status='running', request=config.to_dict(), bases=bases, grids=grids, geometry=geometry,
                  comparison='adjacent input levels, one axis at a time; energy/frequency/distance rank',
                  runs=[], differences=[])
    successes = {}
    save_tables(directory, report)
    try:
        for (basis, grid), cfg in configs.items():
            if control and control.requested:
                report['status'] = 'cancelled'
                break
            progress(f'{cfg.formula} · {cfg.calculation} · {basis} · grid {grid}')
            started = perf_counter()
            row = dict(basis=basis, grid=grid)
            try:
                summary, path = run_calculation(cfg, directory/'runs', backend, control=control)
                row.update(status=summary['status'], directory=str(Path(path).relative_to(directory)))
            except CalculationError as exc:
                row.update(status='error', directory=str(Path(exc.directory).relative_to(directory)), error=str(exc))
            row['wall_seconds'] = perf_counter()-started
            report['runs'].append(row)
            if row['status'] == 'converged':
                values = measure(summary, cfg.calculation)
                current = dict(basis=basis, grid=grid, values=values)
                row['observables'] = [dict(observable=k[0], rank=k[1], unit=k[2], value=v) for k,v in values.items()]
                successes[basis,grid] = current
                i, j = bases.index(basis), grids.index(grid)
                for axis, neighbor in [('basis', (bases[i-1], grid) if i else None),
                                       ('grid', (basis, grids[j-1]) if j else None)]:
                    if neighbor in successes:
                        report['differences'].extend(difference(axis, successes[neighbor], current))
            save_tables(directory, report)
            if row['status'] == 'cancelled' or (control and control.requested):
                report['status'] = 'cancelled'
                break
        else:
            report['status'] = 'completed' if len(successes) == len(configs) else 'incomplete'
    except Exception as exc:
        report.update(status='error', error=str(exc))
        raise
    finally:
        save_tables(directory, report)
    return report, directory


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--config', type=Path)
    source.add_argument('--xyz', type=Path, help='중성 구조 Å. 전하·스핀 등 상세 설정에는 --config 사용')
    parser.add_argument('--calculation', choices=OPERATIONS)
    parser.add_argument('--scan', help='spin_scan 거리 설정')
    parser.add_argument('--bases', nargs='+', default=DEFAULT_BASES)
    parser.add_argument('--grids', nargs='+', type=int, default=DEFAULT_GRIDS)
    parser.add_argument('--backend', choices=('auto','local','wsl'), default='auto')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    try:
        if args.config:
            cfg = load_dft_config(args.config)
            changes = {}
            if args.calculation:
                changes['calculation'] = args.calculation
                if args.calculation != 'spin_scan':
                    changes['scan_spec'] = ''
            if args.scan:
                changes['scan_spec'] = args.scan
            cfg = replace(cfg, **changes)
        else:
            cfg = DFTConfig(atoms=read_xyz(args.xyz), engine='native', calculation=args.calculation or 'uvvis',
                            scan_spec=args.scan or '', cube_points=12, threads=1)
        report, directory = run_comparison(cfg, args.bases, args.grids, args.output, args.backend)
    except (ValueError, RuntimeError, OSError) as exc:
        parser.exit(1, f'비교 오류: {exc}\n')
    print(f"{report['status']} · {directory}")
    return 0 if report['status'] == 'completed' else 2


if __name__ == '__main__':
    raise SystemExit(main())
