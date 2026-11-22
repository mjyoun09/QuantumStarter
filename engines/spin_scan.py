"""Frozen-geometry bond scans and real external spin stability of native closed-shell KS."""
import csv
from dataclasses import replace
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.linalg import eigh

from molecular import Atom, parse_scan, save_dft_config
from scf_control import check_cancelled
from .native_scf import NativeKS
from .native_stability import OrbitalRotations, lowest_mode
from .pyscf_engine import BOHR_ANGSTROM, build_molecule, calculate_molecule, configure_grids, publish_results


class SpinRotations:
    """P.T H P, where P maps v to opposite unit-normalized alpha/beta rotations."""
    def __init__(self, restricted):
        cfg = replace(restricted.config, method='uks', scf_solver='newton', native_stability='none')
        mf = NativeKS(restricted.mol, cfg, restricted.directory)
        mf.grids = restricted.grids
        mf.mo_coeff = np.stack((restricted.mo_coeff, restricted.mo_coeff))
        self.full = OrbitalRotations(mf)
        self.dimension = self.full.sizes[0]

    def hessian_vector(self, v):
        a, b = np.split(self.full.hessian_vector(self.expand(v)), 2)
        return (a-b)/np.sqrt(2)

    @staticmethod
    def expand(v):
        return np.concatenate((v, -v))/np.sqrt(2)


def restricted_point(config, directory):
    mf = NativeKS(build_molecule(config), config, directory)
    configure_grids(mf.grids, config)
    mf.kernel()
    if not mf.converged:
        raise RuntimeError('거리 스캔의 RKS가 수렴하지 않았습니다.')
    rotations = SpinRotations(mf)
    if rotations.dimension == 0:
        raise ValueError('스핀 안정성을 검사하려면 가상 궤도가 필요합니다. 더 큰 기저를 선택하세요.')
    value, mode = lowest_mode(rotations, tolerance=1e-9)
    return mf, rotations, float(value), mode


def continued_density(previous, mol):
    from pyscf.scf.addons import project_mo_nr2nr
    oldmol, occupied = previous
    overlap = mol.intor_symmetric('int1e_ovlp')
    channels = []
    for c in occupied:
        p = project_mo_nr2nr(oldmol, c, mol)
        values, vectors = eigh(p.T @ overlap @ p)
        if values.min() <= 1e-10:
            return None
        orth = p @ (vectors / np.sqrt(values))
        channels.append(orth @ orth.T)
    return np.array(channels)


def follow_point(config, directory, previous):
    rks, rotations, value, mode = restricted_point(config, directory)
    seeds = []
    if value < -1e-7:
        angle = min(1., max(.12, np.sqrt(-value)))
        seeds.append(rotations.full.density(rotations.full.rotated(angle*rotations.expand(mode))))
    if previous is not None:
        guess = continued_density(previous, rks.mol)
        if guess is not None:
            seeds.append(guess)
    uks = rotations.full.mf
    best = rks.e_tot
    density = rks.make_rdm1()
    selected, s2, continuation = config, 0., None
    for seed in seeds:
        trial = uks.kernel(dm0=seed)
        if uks.converged and trial < best:
            trial_s2 = float(uks.spin_square()[0])
            if trial_s2 <= 1e-8:
                # Continuation can return to the restricted solution within SCF roundoff.
                continue
            best, density, selected = trial, uks.make_rdm1(), uks.config
            s2 = trial_s2
            continuation = (uks.mol, [c[:, occ > 0].copy() for c, occ in zip(uks.mo_coeff, uks.mo_occ)])
    if value < -1e-7 and selected is config:
        raise RuntimeError(f'스핀 곡률 {value:.5g} Ha가 음수이나 더 낮은 수렴 UKS 해를 찾지 못했습니다.')
    correction = 0.
    if config.dispersion == 'd4':
        from .dispersion import d4_energy
        correction = d4_energy(rks.mol.atom_charges(), rks.mol.atom_coords(), config.charge, config.xc)['energy_hartree']
    row = dict(rks_hartree=rks.e_tot+correction, lowest_found_hartree=best+correction,
               delta_hartree=best-rks.e_tot, spin_curvature_hartree=value, spin_square=s2,
               selected=selected.resolved_method)
    return row, selected, density, continuation


def scan_bond(config, directory, progress=print):
    from pyscf import lib
    started = perf_counter()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    if (directory/'restart_input.npz').exists():
        raise ValueError('거리 스캔은 단일점 체크포인트로 재시작할 수 없습니다.')
    lib.num_threads(config.threads)
    save_dft_config(config, directory/'settings.json')
    _, atoms, values = parse_scan(config.scan_spec, len(config.atoms))
    i, j = (a-1 for a in atoms)
    original = np.array([a.position for a in config.atoms]) * (BOHR_ANGSTROM if config.unit == 'bohr' else 1.)
    axis = original[j]-original[i]
    norm = np.linalg.norm(axis)
    if norm < 1e-9:
        raise ValueError('선택한 두 원자의 초기 위치가 같습니다.')

    def geometry(distance):
        coords = original.copy()
        coords[j] = coords[i] + distance * (axis/norm)
        return replace(config, atoms=tuple(Atom(a.element, tuple(x)) for a, x in zip(config.atoms, coords)),
                       unit='angstrom', calculation='single_point', scan_spec='')

    rows, brackets, continuation = [], [], None
    with (directory/'spin_scan.csv').open('w', newline='', encoding='utf-8') as stream:
        for distance in values:
            check_cancelled(directory)
            row, final_cfg, density, continuation = follow_point(geometry(distance), directory, continuation)
            row = dict(distance_angstrom=distance, **row)
            if not rows:
                writer = csv.DictWriter(stream, fieldnames=row, lineterminator='\n')
                writer.writeheader()
            rows.append(row)
            writer.writerow(row)
            stream.flush()
            progress(f"결합 {distance:.5f} Å · 곡률 {row['spin_curvature_hartree']:+.6g} Ha · {row['selected']}")
    for a, b in zip(rows, rows[1:]):
        lo, hi = sorted((a,b), key=lambda p: p['distance_angstrom'])
        left, right = lo['distance_angstrom'], hi['distance_angstrom']
        sign = np.sign(lo['spin_curvature_hartree'])
        if sign == np.sign(hi['spin_curvature_hartree']):
            continue
        while right-left > .0004:
            mid = (left+right)/2
            _, _, value, _ = restricted_point(geometry(mid), directory)
            if np.sign(value) == sign:
                left = mid
            else:
                right = mid
        brackets.append([left, right])
    summary = calculate_molecule(final_cfg, directory, progress, guess_density=density, publish=False)
    if not summary['converged'] or abs(summary['total_energy_hartree']-rows[-1]['lowest_found_hartree']) > 1e-7:
        raise RuntimeError('마지막 거리점의 선택된 전자 상태를 재현하지 못했습니다.')
    summary['spin_scan'] = dict(atoms=list(atoms), moved_atom=j+1, points=rows, brackets_angstrom=brackets,
                                table_file='spin_scan.csv', energy_geometry='last requested distance')
    summary['calculation'], summary['config'] = 'spin_scan', config.to_dict()
    summary['final_geometry_config'] = final_cfg.to_dict()
    summary['elapsed_seconds'] = perf_counter() - started
    summary['assumptions'][1] = 'frozen-geometry bond scan; final energy and density belong to the last distance'
    publish_results(summary, directory)
    return summary
