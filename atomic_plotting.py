"""He LDA 방사 격자 결과 그림: 밀도, 퍼텐셜 성분, SCF 수렴 이력."""
import numpy as np
from matplotlib.figure import Figure

from atomic_dft import AtomicDFTResult


def create_atomic_figure(result: AtomicDFTResult) -> Figure:
    fig = Figure(figsize=(12, 3.8), dpi=100, constrained_layout=True, facecolor='white')
    density_ax, potential_ax, scf_ax = fig.subplots(1, 3)
    r = result.r_bohr
    view = r <= 5.0

    density_ax.plot(r[view], result.radial_density[view], color='#2563eb', lw=1.8, label='4πr² n(r)')
    density_ax.plot(r[view], result.u[view], color='#0d9488', lw=1.2, ls='--', label='u₁ₛ(r)')
    density_ax.set_title('Electron density (2 electrons in 1s)', fontsize=10)
    density_ax.set_ylabel('electrons / bohr   ·   bohr^(-1/2)', fontsize=9)

    p = result.potentials
    for key, label, color in (('v_ext', 'v_ext = -2/r', '#64748b'), ('v_hartree', 'v_H', '#d97706'),
                              ('v_xc', 'v_xc (LDA)', '#7c3aed'), ('v_eff', 'v_eff', '#e11d48')):
        potential_ax.plot(r[view], p[key][view], color=color, lw=1.5, label=label)
    potential_ax.axhline(result.eigenvalue_1s, color='#2563eb', lw=0.9, ls=':',
                         label=f'ε₁ₛ = {result.eigenvalue_1s:.6f} Ha')
    potential_ax.set_ylim(-6.0, 3.0)
    potential_ax.set_title('Kohn-Sham potential (last input density)', fontsize=10)
    potential_ax.set_ylabel('potential [hartree]', fontsize=9)
    for ax in (density_ax, potential_ax):
        ax.set_xlabel('r [bohr]', fontsize=9)
        ax.axhline(0, color='#94a3b8', lw=0.6)

    iterations = [h['iteration'] for h in result.history]
    residual = [h['density_residual_electrons'] for h in result.history]
    delta = [abs(h['delta_energy_hartree']) if h['delta_energy_hartree'] else np.nan for h in result.history]
    scf_ax.semilogy(iterations, residual, 'o-', ms=3, color='#2563eb', label='density residual [e]')
    scf_ax.semilogy(iterations, delta, 's-', ms=3, color='#e11d48', label='|ΔE| [Ha]')
    scf_ax.axhline(result.config.density_tol_electrons, color='#2563eb', lw=0.8, ls=':')
    scf_ax.axhline(result.config.energy_tol_hartree, color='#e11d48', lw=0.8, ls=':')
    scf_ax.set_xlabel('SCF iteration', fontsize=9)
    scf_ax.set_title(f'SCF: {result.status} after {len(result.history)} iterations', fontsize=10)
    for ax in (density_ax, potential_ax, scf_ax):
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, fontsize=8)
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle(f'He Kohn-Sham LDA (Slater + VWN5) | E_total = {result.total_energy:.10f} Ha | '
                 f'dr = {result.dr:.3g} bohr, r_max = {result.config.r_max_bohr:g} bohr',
                 fontsize=10, color='#475569')
    return fig
