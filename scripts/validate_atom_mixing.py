"""Pair linear and Anderson radial PBE runs on one thread, writing one JSON line per pair."""
import argparse
from collections import Counter
from dataclasses import asdict, replace
import json
import os
from pathlib import Path
import platform
import sys
from time import perf_counter

if __name__ == '__main__':
    for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
        os.environ[variable] = '1'

import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from radial_atom import AtomConfig, SYMBOLS, quality_warnings, run_atom


PILOT = (('He', False), ('B', True), ('Fe', True), ('Kr', False))


def compare(config, methods):
    runs, metrics = {}, {}
    for method in methods:
        selected = replace(config, mixing_method=method)
        start = perf_counter()
        result = run_atom(selected)
        seconds = perf_counter() - start
        runs[method] = result
        last = result.history[-1]
        metrics[method] = dict(iterations=len(result.history), seconds=seconds, converged=result.converged,
                               energy_hartree=result.total_energy, residual_electrons=last['density_residual_electrons'],
                               delta_hartree=last['delta_energy_hartree'], electron_count=result.electron_count,
                               warnings=quality_warnings(result))
    linear, anderson = runs['linear'], runs['anderson']
    energy = abs(anderson.total_energy - linear.total_energy)
    density = float(np.sum(np.abs(anderson.spin_radial - linear.spin_radial) * linear.grid.weights))
    eigenvalue = max(abs(a['eigenvalue_hartree'] - b['eigenvalue_hartree'])
                     for a, b in zip(linear.orbitals, anderson.orbitals) if a['occupation'] > 0)
    first_identical = all(anderson.history[0][k] == value for k, value in linear.history[0].items())
    return dict(element=config.element, polarized=config.spin_polarized, order=methods, **metrics,
                first_identical=first_identical, abs_energy_hartree=energy, spin_density_l1_electrons=density,
                max_occupied_eigenvalue_difference_hartree=eigenvalue,
                steps=Counter(h['input_mixing_step'] for h in anderson.history),
                passed=bool(linear.converged and anderson.converged and first_identical and
                            energy < 1e-9 and density < 2 * config.density_tol_electrons and
                            abs(anderson.electron_count - config.nuclear_charge) < 1e-9 and
                            not any(m['warnings'] for m in metrics.values())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='New JSONL file')
    parser.add_argument('--representative', action='store_true', help='Only He, B, Fe and Kr')
    parser.add_argument('--repeats', type=int, default=1)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    cases = PILOT if args.representative else [(s, spin) for s in SYMBOLS for spin in (False, True)]
    defaults = AtomConfig(functional='pbe')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    failures = 0
    with args.out.open('x', encoding='utf-8') as stream:
        metadata = dict(platform=platform.platform(), processor=platform.processor(), python=platform.python_version(),
                        numpy=np.__version__, scipy=scipy.__version__, threads=1, config=asdict(defaults),
                        cases=len(cases), repeats=args.repeats)
        print(json.dumps(dict(metadata=metadata)), file=stream, flush=True)
        # Exclude first-use imports and library initialization from the timings.
        for method in ('linear', 'anderson'):
            run_atom(replace(defaults, mixing_method=method, max_iterations=3))
        for repeat in range(args.repeats):
            for index, (element, spin) in enumerate(cases):
                order = ('anderson', 'linear') if (index + repeat + index // 2) % 2 else ('linear', 'anderson')
                pair = compare(replace(defaults, element=element, spin_polarized=spin), order)
                print(json.dumps(dict(repeat=repeat + 1, **pair), allow_nan=False), file=stream, flush=True)
                print(f"{repeat + 1} {element} polarized={spin} {pair['linear']['iterations']} -> "
                      f"{pair['anderson']['iterations']} pass={pair['passed']}", flush=True)
                failures += not pair['passed']
    return int(failures > 0)


if __name__ == '__main__':
    raise SystemExit(main())
