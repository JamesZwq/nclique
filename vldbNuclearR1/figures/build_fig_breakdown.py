"""Phase-breakdown figure: 2 rows x 3 cols.

  Row 1 = build + peel time (ms, log y), 4 lines per panel
  Row 2 = peak memory (MB, log y), 2 lines per panel

Each column is one graph (com-youtube, web-Stanford, web-it-2004).
Data hardcoded from the medians reported in tab:breakdown so the figure
is consistent with the prose; CSV bench_breakdown.csv is the source of
truth for these numbers.
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
    fig, axes = plt.subplots(2, len(GRAPHS), figsize=(13, 3.6), sharex=False)

    for col, g in enumerate(GRAPHS):
        rows = DATA[g]
        ss        = [r[0] for r in rows]
        spin_b    = [r[1] for r in rows]
        spin_p    = [r[2] for r in rows]
        spin_mem  = [r[3] for r in rows]
        cnd_b     = [r[4] for r in rows]
        cnd_p     = [r[5] for r in rows]
        cnd_mem   = [r[6] for r in rows]

        # ----- top row: build + peel time -----
        ax = axes[0][col]
        ax.plot(ss, spin_b, color=SPIN_COLOR, marker="o", markersize=4.5,
                linewidth=1.6, linestyle="-",  label=r"SPIN$^{\star}$ build")
        ax.plot(ss, spin_p, color=SPIN_COLOR, marker="o", markersize=4.5,
                linewidth=1.4, linestyle="--", label=r"SPIN$^{\star}$ peel")
        ax.plot(ss, cnd_b,  color=CND_COLOR,  marker="^", markersize=4.5,
                linewidth=1.6, linestyle="-",  label=r"CND build")
        ax.plot(ss, cnd_p,  color=CND_COLOR,  marker="^", markersize=4.5,
                linewidth=1.4, linestyle="--", label=r"CND peel")
        ax.set_yscale("log")
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
        ax.tick_params(axis="both", which="major", labelsize=8.5)
        ax.tick_params(axis="both", which="minor", labelsize=0)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.set_title(g, fontsize=10.5, fontweight="bold")
        if col == 0:
            ax.set_ylabel("time (ms)", fontsize=10)

        # ----- bottom row: memory -----
        ax = axes[1][col]
        ax.plot(ss, spin_mem, color=SPIN_COLOR, marker="o", markersize=4.5,
                linewidth=1.6, label=r"SPIN$^{\star}$")
        ax.plot(ss, cnd_mem,  color=CND_COLOR,  marker="^", markersize=4.5,
                linewidth=1.6, label="CND")
        ax.set_yscale("log")
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
        ax.tick_params(axis="both", which="major", labelsize=8.5)
        ax.tick_params(axis="both", which="minor", labelsize=0)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.set_xlabel(r"$\boldsymbol{s}$", fontsize=10.5, fontweight="bold")
        if col == 0:
            ax.set_ylabel("memory (MB)", fontsize=10)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=9.5,
               frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot(OUT / "fig_breakdown.pdf")
