"""Case-study cross-r figure: 2 panels (quality, time) x bars grouped by r.

Each (graph, s) cell on the x-axis shows three bars colored by
r in {1, 2, 3}.  The quality panel (left, linear) reads as roughly
equal bars across r within each cell - r=1 matches r=2 and r=3.
The time panel (right, log) shows r=1 dramatically shorter than
r=3 within each cell - r=1 is 67-178x faster.

Data hardcoded from tab:cs4 (com-dblp s=10) and tab:cs7 (com-youtube
s in {5, 10, 15}).
"""
from pathlib import Path

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).parent

# (label, rec50_r1, rec50_r2, rec50_r3, time_r1, time_r2, time_r3)
DATA = [
    ("dblp, s=10",     114, 113, 114,   13.4,  123.1, 2379.6),
    ("youtube, s=5",   283, 286, 288,  274.0, 2580.0, 13112.0),
    ("youtube, s=10",   49,  49,  49,   31.0,  319.0,  2093.0),
    ("youtube, s=15",    0,   0,   0,    9.0,   55.0,    87.0),
]

R_COLORS = {1: "#1f78b4", 2: "#ff7f00", 3: "#e31a1c"}


def plot(out_pdf):
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 2.7))

    labels = [r[0] for r in DATA]
    xs = np.arange(len(DATA))
    width = 0.27

    # ---- panel 1: quality (rec50) ----
    ax = axes[0]
    for ri, r in enumerate([1, 2, 3]):
        ys = [row[1 + ri] for row in DATA]
        ax.bar(xs + (ri - 1) * width, ys, width,
               color=R_COLORS[r], edgecolor="white", linewidth=0.4,
               label=f"$r{{=}}{r}$")
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel(r"rec50 (\# of communities)", fontsize=10)
    ax.set_title("ranking quality (higher is better)",
                 fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=8.5)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.grid(True, axis="y", which="major", alpha=0.25, linestyle=":")
    ax.legend(fontsize=9, frameon=False, loc="upper right", ncol=3)

    # ---- panel 2: time ----
    ax = axes[1]
    for ri, r in enumerate([1, 2, 3]):
        ys = [row[4 + ri] for row in DATA]
        ax.bar(xs + (ri - 1) * width, ys, width,
               color=R_COLORS[r], edgecolor="white", linewidth=0.4,
               label=f"$r{{=}}{r}$")
    ax.set_yscale("log")
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("decomposition time (ms)", fontsize=10)
    ax.set_title("decomposition time (lower is better)",
                 fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=8.5)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.grid(True, axis="y", which="major", alpha=0.25, linestyle=":")
    ax.legend(fontsize=9, frameon=False, loc="upper right", ncol=3)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot(OUT / "fig_crossr.pdf")
