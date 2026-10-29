"""화면에 넣는 그래프와 PNG 출력에 동일한 Figure를 사용한다."""
from pathlib import Path
import numpy as np
from matplotlib.figure import Figure
from solver import BOHR_TO_ANGSTROM, Result


def create_figure(result: Result, selected_state=None, unit="bohr", compare=True) -> Figure:
    fig = Figure(figsize=(10, 3.6), dpi=100, constrained_layout=True, facecolor="white")
    axes = fig.subplots(1, 2)
    z, ell = result.config.nuclear_charge, result.config.angular_momentum_l
    factor = BOHR_TO_ANGSTROM if unit == "angstrom" else 1.0
    unit_label = "angstrom" if unit == "angstrom" else "bohr"
    r = result.r_bohr * factor
    indices = range(result.config.num_states) if selected_state is None else [selected_state]
    for j in indices:
        n = result.principal_n[j]
        color = ("#2563eb", "#0d9488", "#d97706", "#7c3aed", "#e11d48", "#64748b")[j]
        label = f"n={n}"
        axes[0].plot(r, result.u[:, j] / np.sqrt(factor), color=color, lw=1.8, label=label)
        axes[1].plot(r, result.radial_probability[:, j] / factor, color=color, lw=1.8, label=label)
        if compare:
            axes[0].plot(r, result.u_exact[:, j] / np.sqrt(factor), "--", color=color, alpha=0.65, lw=1)
            axes[1].plot(r, np.abs(result.u_exact[:, j]) ** 2 / factor, "--", color=color, alpha=0.65, lw=1)
    axes[0].set_title("Reduced radial wavefunction", fontsize=10, pad=10)
    axes[0].set_ylabel(f"u(r) [{unit_label}^(-1/2)]", fontsize=9)
    axes[0].axhline(0, color="#94a3b8", lw=0.7)
    axes[1].set_title("Radial probability density", fontsize=10, pad=10)
    axes[1].set_ylabel(f"P(r) [{unit_label}^(-1)]", fontsize=9)
    view_max = min(result.config.r_max_bohr, 4.0 * int(result.principal_n[-1]) ** 2 / z) * factor
    for ax in axes:
        ax.set_xlim(0, view_max)
        ax.set_xlabel(f"Distance r [{unit_label}]", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"One electron | Z={z}, l={ell} | solid: numerical"
                 + (", dashed: exact" if compare else ""), fontsize=10, color="#475569")
    return fig


def save_plot(result: Result, destination: Path, show=True):
    figure = create_figure(result)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=170)
    if show:
        import tkinter as tk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        root = tk.Tk()
        root.title("QuantumStarter | 계산 결과")
        canvas = FigureCanvasTkAgg(figure, master=root)
        canvas.draw()
        NavigationToolbar2Tk(canvas, root).update()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        root.mainloop()
