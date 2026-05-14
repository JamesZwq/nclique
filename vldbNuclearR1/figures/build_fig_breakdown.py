"""Phase-breakdown figure: 1 row x 3 cols.

Each column is one graph (com-youtube, web-Stanford, web-it-2004).
Each panel shows two SPIN-star lines on a log-time axis:
  - build (solid)  -- CPI construction time
  - peel (dashed)  -- post-build peel time

The story is the bottleneck-shift inside SPIN-star itself: build stays
near-flat across s while peel drops by 1-3 orders of magnitude.
CND comparison is reported separately in fig:exp-time / fig:hier;
memory in fig:exp-mem.  Data hardcoded from tab:breakdown medians.
"""
from pathlib import Path

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt

OUT = Path(__file__).parent

# Per-graph medians (matched cells from tab:breakdown).
# Schema: graph -> [(s, spin_build, spin_peel, spin_mem,
#                       cnd_build, cnd_peel, cnd_mem), ...]
DATA = {
    "com-youtube": [
        ( 3, 3863,   440.1, 381, 3982,   789.7, 729),
        ( 5, 3730,   273.5, 315, 3795,   581.6, 642),
        ( 8, 2991,    86.3, 187, 3053,   270.2, 460),
        (12, 2701,    13.0, 119, 2649,   148.5, 385),
        (16, 2458,     8.2, 114, 2488,   131.8, 380),
    ],
    "web-Stanford": [
        ( 3, 6602,   381.6, 338, 6548,   521.5, 508),
        (10, 5840,   130.4, 176, 5961,   210.1, 261),
        (20, 5482,    28.9,  97, 5591,    77.8, 151),
        (40, 5087,     9.6,  70, 5007,    49.6, 129),
        (60, 4900,     1.4,  53, 4907,    31.3, 111),
    ],
    "web-it-2004": [
        (  3, 124074, 351.0, 425, 125347, 465.9, 556),
        ( 30, 122558,  63.9, 230, 122689, 173.7, 314),
        (100, 121891,  40.1, 207, 121372, 135.9, 296),
        (200, 113825,  26.1, 178, 107516, 108.5, 273),
        (400,   7850,   4.6, 139,   7847,  62.7, 226),
    ],
}

GRAPHS = ["com-youtube", "web-Stanford", "web-it-2004"]

SPIN_COLOR = "#1f78b4"
CND_COLOR  = "#e31a1c"


def plot(out_pdf):
    fig, axes = plt.subplots(1, len(GRAPHS), figsize=(13, 2.4), sharex=False)

    for col, g in enumerate(GRAPHS):
        rows = DATA[g]
        ss     = [r[0] for r in rows]
        spin_b = [r[1] for r in rows]
        spin_p = [r[2] for r in rows]

        ax = axes[col]
        ax.plot(ss, spin_b, color=SPIN_COLOR, marker="o", markersize=4.5,
                linewidth=1.6, linestyle="-",  label="build")
        ax.plot(ss, spin_p, color=SPIN_COLOR, marker="o", markersize=4.5,
                linewidth=1.4, linestyle="--", label="peel")
        ax.set_yscale("log")
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
        ax.tick_params(axis="both", which="major", labelsize=8.5)
        ax.tick_params(axis="both", which="minor", labelsize=0)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.set_title(g, fontsize=10.5, fontweight="bold")
        ax.set_xlabel(r"$\boldsymbol{s}$", fontsize=10.5, fontweight="bold")
        if col == 0:
            ax.set_ylabel(r"SPIN$^{\star}$ time (ms)", fontsize=10)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, fontsize=9.5,
               frameon=False, bbox_to_anchor=(0.5, 1.04))
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot(OUT / "fig_breakdown.pdf")
