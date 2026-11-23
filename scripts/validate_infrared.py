"""Compare native Cartesian dipole derivatives with PySCF dipoles displaced along normal modes."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

from ase import units
import numpy as np
from pyscf import dft, lib
from pyscf.hessian import thermo

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from engines import calculate
from engines.gas_thermo import write_csv
from engines.pyscf_engine import build_molecule, configure_grids
from molecular import DFTConfig, parse_coordinates
from molecular_isotopes import WATER_EXAMPLE, compare_isotopes
from workflow import new_run_dir


def direct_mode_reference(cfg, flat_hessian, masses, step=0.0025):
    """Step Q in bohr sqrt(u). Solve two external SCFs per independently projected mode."""
    mol = build_molecule(cfg)
    hess = flat_hessian.reshape(mol.natm, 3, mol.natm, 3).transpose(0, 2, 1, 3)
    modes = thermo.harmonic_analysis(mol, hess, mass=masses, imaginary_freq=False)
    intensities = []
    for mode in modes['norm_mode']:
        dipoles = []
        for sign in (1, -1):
            displaced = mol.set_geom_(mol.atom_coords() + sign * step * mode, unit='bohr', inplace=False)
            mf = dft.RKS(displaced)
            mf.xc, mf.small_rho_cutoff = cfg.xc, 0
            configure_grids(mf.grids, cfg)
            mf.conv_tol, mf.conv_tol_grad = 1e-12, 1e-9
            mf.kernel()
            if not mf.converged:
                raise RuntimeError('독립 PySCF 모드 변위 계산이 수렴하지 않았습니다.')
            dipoles.append(mf.dip_moment(unit='AU', origin=(0, 0, 0), verbose=0))
        derivative = (dipoles[0] - dipoles[1]) / (2 * step)
        # ASE uses 42.255 km/mol per (Debye/Å)²/u. Its Debye is expressed in e Å.
        intensities.append(float(np.dot(derivative, derivative) * 42.255 / units.Debye**2))
    return modes['freq_wavenumber'], intensities


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/ir_validation')
    args = parser.parse_args()
    out = new_run_dir(args.output).resolve()
    print(out, flush=True)
    lib.num_threads(1)
    settings = json.loads(WATER_EXAMPLE.read_bytes())['config']
    cfg = DFTConfig.from_dict(dict(settings, engine='native', calculation='frequencies', cube_points=12,
                                   threads=1, conv_tol=1e-11, native_density_tol=1e-8, native_residual_tol=1e-8))
    water = calculate(cfg, out/'water', lambda _: None)
    print('Water, h = 0.005 bohr', flush=True)
    half = calculate(replace(cfg, freq_step_bohr=cfg.freq_step_bohr / 2), out/'water_halfstep', lambda _: None)
    h2 = calculate(replace(cfg, atoms=parse_coordinates('H 0 0 -.375; H 0 0 .375'), dispersion='none'),
                   out/'h2', lambda _: None)
    block = water['frequencies']
    report = compare_isotopes(water, deuterate=True)
    other = compare_isotopes(half, deuterate=True)
    data = []
    for label, molecule in [('reference', 'H2O'), ('substituted', 'D2O')]:
        masses = np.array([a[f'{label}_mass_amu'] for a in report['atoms']])
        freqs, ir = direct_mode_reference(cfg, np.array(block['hessian_hartree_per_bohr2']), masses)
        by_mode = sorted(report['modes'], key=lambda r: r[f'{label}_mode'])
        half_modes = sorted(other['modes'], key=lambda r: r[f'{label}_mode'])
        for i, (row, half_row) in enumerate(zip(by_mode, half_modes)):
            data.append(dict(isotopologue=molecule, mode=i + 1, frequency_cm1=row[f'{label}_cm1'],
                             intensity_km_mol=row[f'{label}_ir_km_mol'], halfstep_cm1=half_row[f'{label}_cm1'],
                             halfstep_km_mol=half_row[f'{label}_ir_km_mol'], pyscf_cm1=float(freqs[i]),
                             pyscf_directional_km_mol=ir[i]))
    write_csv(out/'water_ir_comparison.csv', data)
    # Compact, self-contained desktop example without cube files or displacements.
    example = {k: water[k] for k in ('program', 'version', 'engine', 'system', 'formula', 'calculation',
                                       'converged', 'config', 'warnings')}
    example['frequencies'] = {k: block[k] for k in ('method', 'step_bohr', 'masses_amu', 'frequencies_cm1',
        'reference_max_gradient_hartree_per_bohr', 'hessian_hartree_per_bohr2', 'infrared')}
    (out/'water_ir_reference.json').write_text(json.dumps(example, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    for row in data:
        print(row, flush=True)
    h2_ir, = h2['frequencies']['infrared']['intensities_km_mol']
    print(f'H2 intensity = {h2_ir:.6g} km/mol', flush=True)
    assert h2_ir < 1e-8 and len(data) == 6 and all(r['intensity_km_mol'] > 0 for r in data)
    for r in data:
        assert abs(r['pyscf_cm1'] - r['frequency_cm1']) < .001
        assert abs(r['pyscf_directional_km_mol'] - r['intensity_km_mol']) < .05
        assert abs(r['halfstep_km_mol'] - r['intensity_km_mol']) < .05
        assert abs(r['halfstep_cm1'] - r['frequency_cm1']) < .5


if __name__ == '__main__':
    main()
