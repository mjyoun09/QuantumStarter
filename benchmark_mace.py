"""Matched checkpoint/float64/GPU benchmark: baseline e3nn versus cuEquivariance."""
import argparse
from pathlib import Path
import time
import numpy as np
from ase.io import read
from surface_workflow import dump_json, make_calculator, structure_hash


def benchmark(paths, output, repeats=5):
    import torch
    if repeats < 3:
        raise ValueError('측정은 최소 3회 반복하세요.')
    structures = [read(str(p), format='extxyz') for p in paths]
    results = {}
    predictions = {}
    for acceleration in ('none', 'cueq'):
        calc, metadata = make_calculator('mace', 'cuda', acceleration=acceleration)
        measurements, energies, forces = [], [], []
        for atoms in structures:
            atoms.calc = calc
            # Warm up this shape; measured evaluations explicitly bypass ASE's cache.
            for _ in range(2):
                calc.reset()
                atoms.get_potential_energy()
                atoms.get_forces()
            samples = []
            for _ in range(repeats):
                calc.reset()
                torch.cuda.synchronize()
                start = time.perf_counter()
                energy = float(atoms.get_potential_energy())
                force = atoms.get_forces()
                torch.cuda.synchronize()
                samples.append(time.perf_counter() - start)
            measurements.append(samples)
            energies.append(energy)
            forces.append(force)
        results[acceleration] = dict(model=metadata, samples_seconds=measurements,
                                    median_seconds=float(np.median(measurements)))
        predictions[acceleration] = (np.asarray(energies), np.asarray(forces))
        del calc
    de = float(np.abs(predictions['none'][0] - predictions['cueq'][0]).max())
    df = float(np.abs(predictions['none'][1] - predictions['cueq'][1]).max())
    same_checkpoint = results['none']['model']['checkpoint_sha256'] == results['cueq']['model']['checkpoint_sha256']
    report = dict(results=results, structures=[dict(path=str(p), sha256=structure_hash(a), n_atoms=len(a)) for p, a in zip(paths, structures)],
                  repeats=repeats, speedup=results['none']['median_seconds'] / results['cueq']['median_seconds'],
                  max_energy_difference_ev=de, max_force_component_difference_ev_a=df,
                  same_checkpoint=same_checkpoint, equivalent_within_tolerance=bool(same_checkpoint and de < 1e-6 and df < 1e-5),
                  scope='Same MACE-MP medium checkpoint, float64, GPU, fixed structures; warm energy+force evaluations only. Not DFT or cross-model accuracy.')
    dump_json(report, output)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('structures', nargs='+', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    result = benchmark(args.structures, args.output, args.repeats)
    print({k: v for k, v in result.items() if k not in ('results', 'structures')}, flush=True)
    raise SystemExit(0 if result['equivalent_within_tolerance'] else 2)
