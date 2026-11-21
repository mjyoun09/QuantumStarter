"""H2 PBE/def2-SVP dissociation and real RKS -> UKS curvature, with PySCF checks."""
import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
import pyscf
from pyscf import dft, gto, lib
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from engines.native_scf import NativeKS
from engines.native_stability import OrbitalRotations
from engines.pyscf_engine import configure_grids
from molecular import Atom, DFTConfig


DISTANCES = (.55, .65, .74, .9, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.60, 1.61, 1.62, 1.63, 1.7, 2., 2.5, 3., 4.)


def solvers(config, out):
    mol = gto.M(atom=[(a.element, a.position) for a in config.atoms], basis=config.basis,
                spin=config.spin, unit='Angstrom', verbose=0)
    native = NativeKS(mol, config, out)
    reference = dft.RKS(mol) if config.resolved_method == 'rks' else dft.UKS(mol)
    reference.xc = config.xc
    reference.conv_tol = config.conv_tol
    reference.conv_tol_grad = config.native_residual_tol
    reference.max_cycle = config.max_cycle
    reference.small_rho_cutoff = 0
    reference.chkfile = None
    for solver in (native, reference):
        configure_grids(solver.grids, config)
    return native, reference


def scf(solver, density=None):
    energy = float(solver.kernel(dm0=density))
    if not solver.converged or not np.isfinite(energy):
        raise RuntimeError(f'{type(solver).__name__} SCF did not converge')
    return energy


def geometry(config, distance):
    return replace(config, atoms=(Atom('H', (0., 0., -distance/2)), Atom('H', (0., 0., distance/2))))


def restricted_solution(config, out):
    restricted, reference = solvers(config, out)
    energy = scf(restricted)
    unrestricted, reference_u = solvers(replace(config, method='uks', scf_solver='newton'), out)
    unrestricted.grids = restricted.grids
    unrestricted.mo_coeff = np.stack([restricted.mo_coeff, restricted.mo_coeff])
    rotations = OrbitalRotations(unrestricted)
    # P maps a unit virtual-occupied rotation to (+v,-v)/sqrt(2) in UKS space.
    dimension = rotations.sizes[0]
    projection = np.vstack((np.eye(dimension), -np.eye(dimension))) / np.sqrt(2.)
    hessian = np.column_stack([projection.T @ rotations.hessian_vector(v) for v in projection.T])
    if np.max(np.abs(hessian-hessian.T)) > 1e-7:
        raise RuntimeError('Native spin-antisymmetric Hessian is not symmetric')
    values, vectors = np.linalg.eigh((hessian + hessian.T) / 2)
    return energy, float(values[0]), projection @ vectors[:, 0], rotations, reference, reference_u


def reference_curvature(restricted, unrestricted):
    """Full spin-antisymmetric Hessian by central differences of PySCF UKS gradients."""
    c = restricted.mo_coeff
    n = c.shape[1] - 1  # one occupied spatial orbital for H2
    step = 1e-4
    columns = []
    for j in range(n):
        generator = np.zeros((n + 1, n + 1))
        generator[j+1, 0], generator[0, j+1] = 1/np.sqrt(2.), -1/np.sqrt(2.)
        gradients = []
        for t in (+step, -step):
            cs = np.array([c @ expm(sign*t*generator) for sign in (1., -1.)])
            dm = cs[:, :, :1] @ cs[:, :, :1].transpose(0, 2, 1)
            fock = unrestricted.get_fock(dm=dm)
            g = [2 * spin_c[:, 1:].T @ f @ spin_c[:, 0] for spin_c, f in zip(cs, fock)]
            gradients.append((g[0]-g[1]) / np.sqrt(2.))
        columns.append((gradients[0]-gradients[1])/(2*step))
    matrix = np.column_stack(columns)
    return float(np.linalg.eigvalsh((matrix + matrix.T)/2)[0])


def scan_point(config, out):
    er, curvature, mode, rotations, ref_r, ref_u = restricted_solution(config, out)
    # Opposite spin rotation gives both branches access to broken symmetry.
    angle = min(1., max(.12, np.sqrt(max(0., -curvature))))
    seed = rotations.density(rotations.rotated(angle*mode))
    native_u = rotations.mf
    eu = scf(native_u, seed)
    ref_er = scf(ref_r)
    # Independent reference seed: one normalized 1s AO on each atom, opposite spins.
    indices = [s[2] for s in ref_u.mol.aoslice_by_atom()]
    atomic = np.eye(ref_u.mol.nao_nr())[:, indices]
    atomic /= np.sqrt(np.diag(atomic.T @ ref_u.get_ovlp() @ atomic))
    ref_dm = np.array([np.outer(a, a) for a in atomic.T])
    ref_eu = scf(ref_u, ref_dm)
    dm = native_u.make_rdm1()
    spin_population = np.diag((dm[0]-dm[1]) @ native_u.overlap)
    spin_h1 = float(sum(spin_population[:ref_u.mol.aoslice_by_atom()[0, 3]]))
    s2 = float(native_u.spin_square()[0])
    ref_s2 = float(ref_u.spin_square()[0])
    if (max(abs(er-ref_er), abs(eu-ref_eu)) > 1e-8 or abs(s2-ref_s2) > 1e-4 or
            (curvature < -1e-4 and eu >= er-1e-10) or (curvature > 1e-4 and abs(er-eu) > 1e-8)):
        raise RuntimeError(f'Branch disagreement: curvature={curvature}, native={(er,eu,s2)}, '
                           f'PySCF={(ref_er,ref_eu,ref_s2)}')
    return dict(rks_hartree=er, uks_hartree=eu, uks_minus_rks_hartree=eu-er, spin_curvature_hartree=curvature,
                abs_mulliken_spin_h1=abs(spin_h1), spin_square=s2, pyscf_rks_hartree=ref_er,
                pyscf_uks_hartree=ref_eu, pyscf_spin_square=ref_s2)


def plot(rows, bracket, path):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    fig = Figure(figsize=(9.3, 3.9), layout='constrained')
    FigureCanvasAgg(fig)
    energy, curvature = fig.subplots(1, 2)
    x = [r['distance_angstrom'] for r in rows]
    for key, label, mark in [('rks_hartree', 'RKS', 'o-'), ('uks_hartree', 'UKS', 'x--')]:
        energy.plot(x, [r[key] for r in rows], mark, ms=4, label=label)
    energy.set_ylabel('Total energy (hartree)')
    energy.legend()
    curvature.plot(x, [r['spin_curvature_hartree'] for r in rows], '.-', color='#126568')
    curvature.axhline(0, color='.4', lw=.8)
    curvature.set_ylabel('Lowest spin curvature (hartree)')
    for ax in (energy, curvature):
        ax.axvline(np.mean(bracket), color='.4', lw=1, ls=':')
        ax.set_xlabel('H–H distance (Å)')
    fig.suptitle('H₂ dissociation   |   PBE / def2-SVP')
    fig.savefig(path, dpi=165)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='new output directory')
    out = parser.parse_args().out
    out.mkdir(parents=True, exist_ok=False)
    lib.num_threads(1)
    config = DFTConfig(engine='native', method='rks', xc='pbe', basis='def2-svp', grid_level=3,
                       max_cycle=350, conv_tol=1e-11, native_density_tol=1e-8, native_residual_tol=1e-8, threads=1)
    rows = []
    with (out/'h2_scan.csv').open('x', newline='', encoding='utf-8') as stream:
        for distance in DISTANCES:
            cfg = geometry(config, distance)
            row = dict(distance_angstrom=distance, **scan_point(cfg, out))
            if not rows:
                writer = csv.DictWriter(stream, fieldnames=row.keys(), lineterminator='\n')
                writer.writeheader()
            writer.writerow(row)
            stream.flush()
            rows.append(row)
            print(f"R={distance:g} A, UKS-RKS={row['uks_minus_rks_hartree']:.7f} Ha, "
                  f"spin curvature={row['spin_curvature_hartree']:.6f}", flush=True)
    brackets = [(a['distance_angstrom'], b['distance_angstrom']) for a, b in zip(rows, rows[1:])
                if a['spin_curvature_hartree'] > 0 > b['spin_curvature_hartree']]
    if len(brackets) != 1:
        raise RuntimeError(f'Expected one positive-to-negative curvature crossing: {brackets}')
    left, right = brackets[0]
    trace = []
    while right-left > .0005:
        mid = (left+right)/2
        _, value, *_ = restricted_solution(geometry(config, mid), out)
        trace.append([mid, value])
        if value > 0:
            left = mid
        else:
            right = mid
    endpoints = []
    for distance in (left, right):
        _, value, _, _, ref_r, ref_u = restricted_solution(geometry(config, distance), out)
        scf(ref_r)
        independent = reference_curvature(ref_r, ref_u)
        if value*independent <= 0 or abs(value-independent) > 1e-6:
            raise RuntimeError(f'PySCF bracket check failed at {distance}: {value}, {independent}')
        _, fine_grid, *_ = restricted_solution(geometry(replace(config, grid_level=4), distance), out)
        endpoints.append(dict(distance_angstrom=distance, native_hartree=value, pyscf_fd_hartree=independent,
                              level4_native_hartree=fine_grid))
    # Same XC/basis on the isolated doublet fragments, with no ghost functions.
    one_h = replace(config, atoms=(Atom('H', (0., 0., 0.)),), spin=1, method='uks')
    atom, ref_atom = solvers(one_h, out)
    atomic_energy, reference_atomic_energy = scf(atom), scf(ref_atom)
    if abs(atomic_energy - reference_atomic_energy) > 1e-8:
        raise RuntimeError('Isolated H reference mismatch')
    report = dict(template=geometry(config, .74).to_dict(), native_uks_solver='newton', pyscf_uks_solver='diis',
                  pyscf=pyscf.__version__, numpy=np.__version__,
                  isolated_h_hartree=atomic_energy, pyscf_isolated_h_hartree=reference_atomic_energy,
                  bracket=endpoints, bisection=trace, max_energy_difference_hartree=max(
                      abs(r[m+'_hartree']-r['pyscf_'+m+'_hartree']) for r in rows for m in ('rks','uks')))
    (out/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    plot(rows, (left, right), out/'h2_dissociation.png')
    print('Curvature crossing (angstrom):', left, right, flush=True)


if __name__ == '__main__':
    main()
