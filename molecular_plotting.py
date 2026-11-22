"""Scientific views of saved samples; plot grids are independent of SCF quadrature."""
from pathlib import Path

import numpy as np
from matplotlib.colors import PowerNorm
from matplotlib.figure import Figure
from scipy.constants import physical_constants

from molecular import ELEMENTS

BOHR_TO_ANGSTROM = physical_constants['Bohr radius'][0] / 1e-10
COLORS = {'H': '#c1c9d7', 'C': '#374151', 'N': '#3b82f6', 'O': '#ef4444',
          'F': '#10b981', 'S': '#eab308', 'P': '#f97316', 'Cl': '#22c55e', 'He': '#a78bfa',
          'K': '#8b5cf6', 'Ca': '#65a30d', 'Sc': '#94a3b8', 'Ti': '#64748b',
          'V': '#a1a1aa', 'Cr': '#0d9488', 'Mn': '#c026d3', 'Fe': '#d97706', 'Co': '#2563eb',
          'Ni': '#059669', 'Cu': '#b45309', 'Zn': '#71717a', 'Ga': '#a78bfa', 'Ge': '#78716c',
          'As': '#a855f7', 'Se': '#f59e0b', 'Br': '#991b1b', 'Kr': '#67e8f9',
          'Rb': '#7c3aed', 'Sr': '#84cc16', 'Y': '#06b6d4', 'Zr': '#0891b2'}


def read_cube(path):
    with Path(path).open(encoding='utf-8') as stream:
        stream.readline()
        stream.readline()
        row = stream.readline().split()
        count, origin = int(row[0]), np.array(row[1:4], dtype=float)
        shape, steps = [], []
        for _ in range(3):
            row = stream.readline().split()
            shape.append(int(row[0]))
            steps.append([float(v) for v in row[1:4]])
        if count < 0 or any(n < 1 for n in shape):
            raise ValueError('이 뷰어는 프로그램이 출력한 직교 CUBE 형식을 지원합니다.')
        for _ in range(count):
            stream.readline()
        values = np.fromstring(stream.read(), sep=' ')
    if values.size != np.prod(shape):
        raise ValueError('CUBE 데이터 크기가 격자 수와 일치하지 않습니다.')
    return values.reshape(shape), origin, np.array(steps) * (np.array(shape) - 1)[:, None]


def analysis_figure(summary):
    fig = Figure(figsize=(9, 4.5), layout='constrained')
    if 'uvvis' in summary:
        peaks = summary['uvvis']['peaks']
        ax = fig.subplots()
        x = [p['wavelength_nm'] for p in peaks]
        y = [p['oscillator_strength'] for p in peaks]
        ax.vlines(x, 0, y, color='#176389', lw=2)
        ax.scatter(x, y, s=16, color='#176389')
        ax.set(xlabel='Wavelength (nm)', ylabel='Oscillator strength', title=f"{summary['formula']} | Singlet TDA")
        ax.set_ylim(bottom=0)
    elif 'gas_thermo' in summary:
        rows = summary['gas_thermo']['rows']
        for ax, key, label in zip(fig.subplots(1, 3), ('cp_J_mol_K','delta_h_kJ_mol','s_J_mol_K'),
                                  ('Cp (J/mol/K)', 'H − H298 (kJ/mol)', 'S (J/mol/K)')):
            for prefix, style, name in [('', '-', 'RRHO'), ('nist_', '--', 'NIST')]:
                ax.plot([r['temperature_K'] for r in rows], [r[prefix+key] for r in rows], style, label=name)
            ax.set(xlabel='Temperature (K)', ylabel=label)
        fig.axes[0].legend()
        fig.suptitle(f"{summary['formula']} | Gas at 1 bar")
    else:
        block = summary['spin_scan']
        rows = block['points']
        x = [r['distance_angstrom'] for r in rows]
        en, cv = fig.subplots(1, 2)
        for key, label, marker in [('rks_hartree','RKS','o-'), ('lowest_found_hartree','Lowest found','x--')]:
            en.plot(x, [r[key] for r in rows], marker, label=label)
        en.set_ylabel('Total energy (hartree)')
        en.legend()
        cv.plot(x, [r['spin_curvature_hartree'] for r in rows], '.-')
        cv.axhline(0, color='.4', lw=.8)
        cv.set_ylabel('Lowest spin curvature (hartree)')
        for ax in (en, cv):
            ax.set_xlabel('Bond distance (Å)')
            for interval in block['brackets_angstrom']:
                ax.axvline(np.mean(interval), color='.5', ls=':')
    return fig


def create_molecular_figure(summary, directory, field=None, plane='yz'):
    if field in (None, '분석 그래프') and any(k in summary for k in ('uvvis', 'gas_thermo', 'spin_scan')):
        return analysis_figure(summary)
    field = field or 'density.cube'
    directory = Path(directory)
    with np.load(directory / 'density_grid.npz', allow_pickle=False) as data:
        coords = data['atom_coords_bohr'].copy()
        numbers = data['atom_numbers'].copy()
        density = data['density'].copy()
        origin = data['origin_bohr'].copy()
        axes = data['axes_bohr'].copy()
    if field not in summary['cube_files'] or Path(field).name != field:
        raise ValueError('저장된 밀도 또는 오비탈 파일을 선택하세요.')
    if field != 'density.cube':
        density, origin, axes = read_cube(directory / field)
    figure = Figure(figsize=(10, 4.5), dpi=100, facecolor='white', layout='constrained')
    atoms = figure.add_subplot(121, projection='3d')
    xyz = coords * BOHR_TO_ANGSTROM
    symbols = [ELEMENTS[int(z) - 1] for z in numbers]
    atoms.scatter(*xyz.T, s=[100 if s == 'H' else 180 for s in symbols],
                  c=[COLORS.get(s, '#a78bfa') for s in symbols], edgecolors='white', depthshade=True)
    for point, symbol in zip(xyz, symbols):
        atoms.text(*point, ' ' + symbol, fontsize=9)
    center = xyz.mean(axis=0)
    radius = max(float(np.ptp(xyz, axis=0).max()) / 2 + 0.5, 1.0)
    for setter, c in zip((atoms.set_xlim, atoms.set_ylim, atoms.set_zlim), center):
        setter(c - radius, c + radius)
    atoms.set_box_aspect((1, 1, 1))
    atoms.set(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]', title=f'{summary["formula"]} | input geometry')
    horizontal, vertical, normal = {'xy': (0, 1, 2), 'xz': (0, 2, 1), 'yz': (1, 2, 0)}[plane]
    lengths = np.diag(axes)
    n = density.shape[normal]
    index = int(np.clip(round((coords.mean(axis=0)[normal] - origin[normal]) / lengths[normal] * (n - 1)), 0, n - 1))
    section = np.take(density, index, axis=normal).T
    normal_value = (origin[normal] + lengths[normal] * index / (n - 1)) * BOHR_TO_ANGSTROM
    extent = np.array([origin[horizontal], origin[horizontal] + lengths[horizontal],
                       origin[vertical], origin[vertical] + lengths[vertical]]) * BOHR_TO_ANGSTROM
    ax = figure.add_subplot(122)
    signed = field != 'density.cube'
    if signed:
        maximum = max(float(np.max(np.abs(section))), 1e-12)
        artist = ax.imshow(section, origin='lower', extent=extent, cmap='RdBu_r',
                           vmin=-maximum, vmax=maximum, interpolation='bilinear')
    else:
        artist = ax.imshow(np.maximum(section, 0), origin='lower', extent=extent, cmap='magma',
                           norm=PowerNorm(gamma=0.35), interpolation='bilinear')
    label = 'e / bohr³' if 'density' in field else 'bohr^(-3/2)'
    figure.colorbar(artist, ax=ax, label=label, shrink=0.76)
    model_label = ('ECP / ' if summary.get('ecp_core_electrons', 0) else
                   'X2C / ' if summary.get('relativity') == 'x2c' else '')
    ax.set(xlabel=f'{plane[0]} [Å]', ylabel=f'{plane[1]} [Å]',
           title=f'{model_label}{field.removesuffix(".cube")}\n{"xyz"[normal]} = {normal_value:.3f} Å')
    if not summary['converged']:
        figure.suptitle('NOT CONVERGED — diagnostic data only', color='#b45309')
    return figure
