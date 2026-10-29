"""Persist gradient checks: frozen-grid finite differences, grid-response size, engine agreement.

Run in WSL: .venv-dft/bin/python scripts/validate_gradient.py
No reference values are hardcoded; every number is computed and written to report.json.
"""
from datetime import datetime
import json
from pathlib import Path
import sys
import tempfile
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pyscf import dft, gto, lib
from molecular import DFTConfig, parse_coordinates
from engines.native_scf import NativeKS
from engines.native_gradient import native_gradient
from engines.pyscf_engine import configure_grids

WATER = 'O 0 0 0; H 0 -.757 .587; H 0 .757 .587'
NH2 = 'N 0 0 0; H 0 -.80 .62; H 0.05 .80 .60'
TIGHT = dict(conv_tol=1e-12, native_density_tol=1e-10, native_residual_tol=1e-9, max_cycle=300)


def solve(mol, cfg, directory, grids=None, dm0=None):
    mf = NativeKS(mol, cfg, directory)
    if grids is None:
        configure_grids(mf.grids, cfg)
    else:
        mf.grids = grids
    energy = mf.kernel(dm0=dm0)
    if not mf.converged:
        raise RuntimeError('native SCF did not converge')
    return mf, energy


def finite_difference(mol, cfg, directory, reference, frozen, step):
    numeric = np.zeros((mol.natm, 3))
    for atom in range(mol.natm):
        for axis in range(3):
            energies = []
            for sign in (1, -1):
                coords = mol.atom_coords().copy()
                coords[atom, axis] += sign * step
                moved = mol.set_geom_(coords, unit='Bohr', inplace=False)
                energies.append(solve(moved, cfg, directory, reference.grids if frozen else None,
                                      reference.make_rdm1())[1])
            numeric[atom, axis] = (energies[0] - energies[1]) / (2 * step)
    return numeric


def main():
    lib.num_threads(2)
    output = ROOT / 'outputs' / 'gradient_validation' / datetime.now().strftime('%Y%m%d_%H%M%S')
    output.mkdir(parents=True)
    report = {'created': datetime.now().isoformat(), 'units': 'hartree/bohr', 'frozen_grid': [], 'grid_response': [],
              'engine_agreement': [], 'pyscf_paths': []}
    with tempfile.TemporaryDirectory() as directory:
        for name, atoms, spin, xc in (('H2O', WATER, 0, 'lda,vwn'), ('H2O', WATER, 0, 'pbe'), ('H2O', WATER, 0, 'pbe0'),
                                      ('NH2', NH2, 1, 'pbe'), ('NH2', NH2, 1, 'b3lyp')):
            cfg = DFTConfig(atoms=parse_coordinates(atoms), spin=spin, xc=xc, basis='sto-3g', engine='native',
                            grid_level=1, **TIGHT)
            mol = gto.M(atom=atoms, basis='sto-3g', spin=spin, verbose=0)
            reference, _ = solve(mol, cfg, directory)
            analytic = sum(native_gradient(reference).values())
            numeric = finite_difference(mol, cfg, directory, reference, True, 1e-4)
            report['frozen_grid'].append(dict(system=name, xc=xc, spin=spin, basis='sto-3g', grid_level=1, step_bohr=1e-4,
                                              max_abs_difference=float(np.abs(analytic - numeric).max()),
                                              max_abs_gradient=float(np.abs(analytic).max())))
            print(report['frozen_grid'][-1], flush=True)
        for level in (1, 2, 3, 4):
            cfg = DFTConfig(atoms=parse_coordinates(WATER), xc='pbe', basis='sto-3g', engine='native', grid_level=level, **TIGHT)
            mol = gto.M(atom=WATER, basis='sto-3g', verbose=0)
            reference, _ = solve(mol, cfg, directory)
            analytic = sum(native_gradient(reference).values())
            numeric = finite_difference(mol, cfg, directory, reference, False, 1e-3)
            report['grid_response'].append(dict(system='H2O', xc='pbe', grid_level=level, points=len(reference.grids.weights),
                                                max_abs_difference_moving_grid=float(np.abs(analytic - numeric).max()),
                                                max_abs_net_force=float(np.abs(analytic.sum(axis=0)).max())))
            print(report['grid_response'][-1], flush=True)
        for basis, xc in (('def2-svp', 'pbe'), ('def2-svp', 'b3lyp')):
            cfg = DFTConfig(atoms=parse_coordinates(WATER), xc=xc, basis=basis, engine='native', grid_level=3, **TIGHT)
            mol = gto.M(atom=WATER, basis=basis, verbose=0)
            own, _ = solve(mol, cfg, directory)
            mine = sum(native_gradient(own).values())
            ref = dft.RKS(mol)
            ref.xc, ref.small_rho_cutoff, ref.conv_tol = xc, 0, 1e-12
            configure_grids(ref.grids, cfg)
            ref.kernel(dm0=own.make_rdm1())
            theirs = ref.nuc_grad_method().kernel()
            report['engine_agreement'].append(dict(system='H2O', basis=basis, xc=xc, grid_level=3,
                                                   max_abs_difference=float(np.abs(mine - theirs).max())))
            print(report['engine_agreement'][-1], flush=True)
    # External-engine paths the GUI can request: Newton SCF and scalar X2C.
    for label, atoms, basis, wrap in (('newton', WATER, 'def2-svp', lambda mf: mf.newton()),
                                      ('sfx2c1e', 'Ag 0 0 0; H 0 0 1.62', 'dyall-v2z', lambda mf: mf.sfx2c1e())):
        mol = gto.M(atom=atoms, basis=basis, verbose=0)
        mf = wrap(dft.RKS(mol))
        mf.xc = 'pbe'
        mf.kernel()
        try:
            gradient = mf.nuc_grad_method().kernel()
            report['pyscf_paths'].append(dict(path=label, supported=True, converged=bool(mf.converged),
                                              max_abs_net_force=float(np.abs(gradient.sum(axis=0)).max())))
        except Exception as exc:  # Record, do not hide, unsupported combinations.
            report['pyscf_paths'].append(dict(path=label, supported=False, error=repr(exc)[:300]))
        print(report['pyscf_paths'][-1], flush=True)
    (output / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('saved', output / 'report.json')


if __name__ == '__main__':
    main()
