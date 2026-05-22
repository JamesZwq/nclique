"""Phase-breakdown figure: 1 row x 6 cols of stacked bars.

Each column is one graph (smallest to largest), each bar one s value.
Bar = build (bottom) + peel (top); linear y-axis so the visual ratio
matches the time ratio.  Peel is a thin sliver on top of a tall build
base, making the "construction dominates total time" story obvious.

Data: paper_data/bench_r1_hierarchy.csv, ST_V3 rows (status OK).
CND comparison lives in fig:exp-time / fig:hier; memory in fig:exp-mem.
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt

OUT      = Path(__file__).parent
DATA_CSV = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data/bench_r1_hierarchy.csv")

# Six representative graphs, smallest to largest by edge count:
# 0.93M / 1.05M / 1.99M / 2.99M / 7.18M / 22.3M.
GRAPHS = [
    ("com-amazon.ungraph",      "com-amazon"),
    ("com-dblp",                "com-dblp"),
    ("web-Stanford",            "web-Stanford"),
    ("com-youtube",             "com-youtube"),
    ("web-it-2004",             "web-it-2004"),
    ("soc-pokec-relationships", "soc-pokec"),
]
S_VALUES = [3, 5, 7, 9, 11, 13, 15]

BUILD_COLOR = "#9ecae1"  # light blue (the bulk)
PEEL_COLOR  = "#ef6548"  # orange-red (the sliver on top)


def load():
    data = defaultdict(dict)   # graph -> {s: (build, peel)}
    with open(DATA_CSV) as f:
        for r in csv.DictReader(f):
            if r["algo"] != "ST_V3": continue
            if r["status"] != "OK": continue
            try:
                s = int(r["s"])
                build = float(r["build_ms"])
                peel  = float(r["peel_ms"])
            except (ValueError, TypeError):
                continue
            if build < 0 or peel < 0: continue
            data[r["graph"]][s] = (build, peel)
    return data


def fmt_pct(pct):
    if pct < 0.1:
        return "<0.1%"
    return f"{pct:.1f}%"


def plot(out_pdf):
    data = load()
    n = len(GRAPHS)
    fig, axes = plt.subplots(1, n, figsize=(n * 2.1, 1.7), sharex=False)
    if n == 1: axes = [axes]

    for col, (stem, label) in enumerate(GRAPHS):
        cells = data.get(stem, {})
        ss     = sorted(s for s in S_VALUES if s in cells)
        if not ss:
            axes[col].set_visible(False); continue
        spin_b = [cells[s][0] for s in ss]
        spin_p = [cells[s][1] for s in ss]
        x_pos  = list(range(len(ss)))

        ax = axes[col]
        ax.bar(x_pos, spin_b, color=BUILD_COLOR, edgecolor="white",
               linewidth=0.4, label="build")
        ax.bar(x_pos, spin_p, bottom=spin_b, color=PEEL_COLOR,
               edgecolor="white", linewidth=0.4, label="peel")

        for i, (b, p) in enumerate(zip(spin_b, spin_p)):
            pct = 100.0 * p / (b + p) if (b + p) > 0 else 0.0
            ax.text(x_pos[i], b + p, fmt_pct(pct), ha="center", va="bottom",
                    fontsize=7.0, color="#404040")

        ax.set_xticks(x_pos)
        ax.set_xticklabels([str(s) for s in ss])
        ax.tick_params(axis="both", which="major", labelsize=8.0)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.grid(True, axis="y", which="major", alpha=0.25, linestyle=":")
        top = max(b + p for b, p in zip(spin_b, spin_p))
        ax.set_ylim(0, top * 1.22)
        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.set_xlabel(r"$\boldsymbol{s}$", fontsize=10, fontweight="bold")
        if col == 0:
            ax.set_ylabel(r"SPIN$^{\star}$ time (ms)", fontsize=9.5)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, fontsize=9.5,
               frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot(OUT / "fig_breakdown.pdf")
