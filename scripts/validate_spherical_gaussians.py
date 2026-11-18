"""Independent spherical Gaussian RKS cross-check for closed-shell atoms.

Uses PySCF's Libxc interface and SciPy Galerkin matrices. Shares no orbital,
quadrature, Coulomb or XC implementation with the finite-difference engine.
Run in WSL with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre
from scipy.linalg import eigh, svd
from pyscf.dft import libxc


# Number of occupied radial functions per angular momentum (all m occupied).
ATOMS = {'He': (2, (1,)), 'Ar': (18, (3, 2)), 'Zn': (30, (4, 2, 1)), 'Kr': (36, (4, 3, 1))}


def quadrature(panels, degree=16):
    edges = np.r_[0, np.geomspace(1e-7, 90, panels)]
    nodes, weights = legendre.leggauss(degree)
    half = np.diff(edges)[:, None] / 2
    r = (edges[:-1, None] + half * (1 + nodes)).ravel()
    # Integrate the Legendre interpolation polynomial from -1 to each node.
    vand = legendre.legvander(nodes, degree - 1)
    primitive = legendre.legint(np.eye(degree), lbnd=-1)
    prefix = legendre.legvander(nodes, degree) @ primitive @ np.linalg.inv(vand)
    return r, (half * weights).ravel(), half, prefix, weights


def cumulative(f, half, prefix, weights):
    samples = f.reshape(len(half), -1)
    panel_integrals = half[:, 0] * (samples @ weights)
    before = np.r_[0, np.cumsum(panel_integrals[:-1])]
    return (before[:, None] + half * (samples @ prefix.T)).ravel()


def gaussian_channel(r, w, ell, z, ratio):
    alpha = .008 * ratio ** np.arange(1 + int(np.ceil(np.log(3e8 / .008) / np.log(ratio))))
    radii = r[:, None]
    b = radii ** ell * np.exp(-radii ** 2 * alpha)
    db = (ell / radii - 2 * radii * alpha) * b
    wr = (w * r ** 2)[:, None]
    # Use the positive kinetic + z² overlap metric. Both the Hamiltonian
    # and its shift stay well scaled despite tight core Gaussians.
    metric = np.vstack((z * np.sqrt(wr) * b, np.sqrt(.5 * wr) * db,
                        np.sqrt(.5 * ell * (ell + 1) * w[:, None]) * b))
    norm = np.linalg.norm(metric, axis=0)
    u, singular, _ = svd(metric / norm, full_matrices=False)
    keep = singular > 3e-9
    basis = u[:len(r), keep] / (z * np.sqrt(wr))
    grad = u[len(r):2 * len(r), keep] / np.sqrt(.5 * wr)
    overlap = basis.T @ (wr * basis)
    kinetic = .5 * (grad.T @ (wr * grad) + ell * (ell + 1) * basis.T @ (w[:, None] * basis))
    core = kinetic - z * basis.T @ (wr / radii * basis)
    return basis, grad, core, overlap


def xc(density, derivative):
    zeros = np.zeros_like(density)
    return libxc.eval_xc('PBE', np.array([density, derivative, zeros, zeros]), deriv=1)[:2]


def run(symbol, ratio, panels):
    z, occupied = ATOMS[symbol]
    r, w, half, prefix, weights = quadrature(panels)
    wr = w * r ** 2
    channels = [gaussian_channel(r, w, ell, z, ratio) for ell in range(len(occupied))]
    # A smooth screened Coulomb guess. n and n' are mixed together thereafter.
    exponent = 2 * z ** (1 / 3)
    n = z * exponent ** 3 / (8 * np.pi) * np.exp(-exponent * r)
    dn = -exponent * n
    energy = None
    for iteration in range(1, 201):
        inner = cumulative(4 * np.pi * r ** 2 * n, half, prefix, weights)
        outer_density = 4 * np.pi * r * n
        vh = inner / r + np.dot(w, outer_density) - cumulative(outer_density, half, prefix, weights)
        _, (vrho, vsigma) = xc(n, dn)
        nout, gout = np.zeros_like(n), np.zeros_like(dn)
        one_body = 0.
        for ell, ((b, db, core, overlap), count) in enumerate(zip(channels, occupied)):
            weak = db.T @ ((wr * 2 * vsigma * dn)[:, None] * b)
            fock = core + b.T @ ((wr * (vh + vrho))[:, None] * b) + weak + weak.T
            # Invert a positive shifted Hamiltonian to resolve low eigenvalues
            # in the presence of very tight Gaussians and high virtual energies.
            dim = len(core)
            inv, vectors = eigh(overlap, (fock + fock.T) / 2 + z ** 2 * overlap,
                                subset_by_index=(dim - count, dim - 1))
            coeff = vectors / np.sqrt(inv)
            orbitals, grads = b @ coeff, db @ coeff
            occupation = 2 * (2 * ell + 1)
            nout += occupation / (4 * np.pi) * np.sum(orbitals ** 2, axis=1)
            gout += occupation / (2 * np.pi) * np.sum(orbitals * grads, axis=1)
            one_body += occupation * np.sum(coeff * (core @ coeff))
        residual = 4 * np.pi * np.dot(wr, np.abs(nout - n))
        charge = 4 * np.pi * r ** 2 * nout
        # Continuous spherical Hartree energy = integral Q(r) dQ(r)/r.
        hartree = np.dot(w, charge / r * cumulative(charge, half, prefix, weights))
        new_energy = one_body + hartree + np.dot(w * charge, xc(nout, gout)[0])
        if energy is not None and abs(new_energy - energy) < 2e-10 and residual < 2e-9:
            break
        energy = new_energy
        n, dn = .7 * n + .3 * nout, .7 * dn + .3 * gout
    else:
        raise RuntimeError(f'{symbol} did not converge: density residual {residual:g}')
    return dict(element=symbol, ratio=ratio, panels=panels, total_energy=float(new_energy),
                iterations=iteration, residual=float(residual), electrons=float(w @ charge))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radial-report', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    reference = json.loads(args.radial_report.read_text(encoding='utf-8'))
    args.out.mkdir(parents=True, exist_ok=False)
    report = {'libxc': libxc.__version__, 'method': 'spherical Gaussian Galerkin RKS with Libxc PBE',
              'tolerance_hartree': 1e-6, 'cases': [], 'complete': False}
    def save():
        path = args.out / 'report.json'
        tmp = path.with_suffix('.tmp')
        tmp.write_text(json.dumps(report, indent=2), encoding='utf-8')
        tmp.replace(path)
    for symbol, (z, _) in ATOMS.items():
        radial = next(c for c in reference['cases'] if c['element'] == symbol and not c['polarized'])
        trials = []
        for ratio, panels in ((1.35, 70), (1.25, 70), (1.25, 100)):
            trial = run(symbol, ratio, panels)
            trials.append(trial)
            print(json.dumps(trial), flush=True)
        energy = trials[-1]['total_energy']
        errors = {'radial_minus_gaussian': radial['richardson'] - energy,
                  'basis_change': trials[1]['total_energy'] - trials[0]['total_energy'],
                  'quadrature_change': energy - trials[1]['total_energy']}
        passed = (radial['converged'] and radial['coarse_converged'] and
                  all(abs(e) < 1e-6 for e in errors.values()) and
                  all(abs(t['electrons'] - z) < 1e-9 for t in trials))
        report['cases'].append(dict(element=symbol, runs=trials, passed=passed, **errors))
        save()
    report.update(complete=True, passed=all(c['passed'] for c in report['cases']))
    save()
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
