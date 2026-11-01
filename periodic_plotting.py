"""Cell geometry and SCF history; no molecular orbital plot for a periodic cell."""
import itertools
import numpy as np
from matplotlib.figure import Figure


def create_periodic_figure(summary):
    config = summary.get('final_geometry_config', summary['config'])
    figure = Figure(figsize=(10, 5), layout='constrained')
    figure.get_layout_engine().set(wspace=.18)
    ax = figure.add_subplot(121, projection='3d')
    lattice = np.asarray(config['lattice'])
    for corner in itertools.product((0, 1), repeat=3):
        for d in range(3):
            if corner[d] == 0:
                end = list(corner)
                end[d] = 1
                p = np.asarray([corner, end]) @ lattice
                ax.plot(*p.T, color='#94a3b8', linewidth=.8)
    for e in sorted({a['element'] for a in config['atoms']}):
        xyz = np.asarray([a['position'] for a in config['atoms'] if a['element'] == e])
        ax.scatter(*xyz.T, s=55, label=e)
    ax.set(xlabel='x / Å', ylabel='y / Å', zlabel='z / Å', title=f"{summary['formula']} · {config['boundary']}")
    ax.legend()
    ax.set_box_aspect(np.maximum(np.ptp(np.asarray(list(itertools.product((0, 1), repeat=3))) @ lattice, axis=0), 1))
    trace = figure.add_subplot(122)
    history = summary.get('scf_history', [])
    if history:
        trace.plot([r['cycle'] for r in history], [r['energy_hartree'] for r in history], '.-')
    trace.set(xlabel='SCF cycle', ylabel='E / Ha per cell', title=summary.get('status', 'Input structure'))
    trace.grid(alpha=.2)
    return figure
