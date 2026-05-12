"""
Scalability figure: SPIN★ time / RSS vs s on five graphs, one curve per
vertex-induced subsample ratio. Single row, one column per graph.

Reads paper_data/scalability_{graph}.csv for each graph in GRAPHS.
"""
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import csv
from collections import defaultdict
from pathlib import Path

OUT      = Path(__file__).parent
DATA_DIR = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data")

# Graphs to include, in size order (small → large). Each entry is
# (csv-stem, display-label).
GRAPHS = [
    ("com-dblp",                 "com-dblp"),
    ("com-youtube",              "com-youtube"),
    ("web-Stanford",             "web-Stanford"),
    ("web-Google",               "web-Google"),
    ("soc-pokec-relationships",  "soc-pokec"),
]

# Only SPIN★ — drop CND so the figure focuses on our algorithm.
ALGO_KEY = "Ours_ST"

# 5 ratios -> 5 colours from a perceptually ordered palette.
RATIO_COLORS = {
    0.2: "#fee08b",
    0.4: "#fdae61",
    0.6: "#f46d43",
    0.8: "#d73027",
    1.0: "#a50026",
}


def load(csv_path):
    """Returns {(algo, ratio): [(s, time_ms, mem_kb), ...] sorted}."""
    if not csv_path.exists():
        return None
    rows = list(csv.DictReader(open(csv_path)))
    data = defaultdict(list)
    for r in rows:
        if r.get("status") != "OK": continue
        try:
            ratio = float(r["ratio"])
            s     = int(r["s"])
            t     = float(r["time_ms"])
            m     = float(r["memory_kB"])
        except (KeyError, ValueError):
            continue
        data[(r["algorithm"], ratio)].append((s, t, m))
    for k in data: data[k].sort()
    return data


def plot(metric, ylabel, out_pdf, mem=False):
    panels = []
    for stem, label in GRAPHS:
        d = load(DATA_DIR / f"scalability_{stem}.csv")
        if d:
            panels.append((label, d))
    if not panels: return

    n_cols = len(panels)
    fig, axes = plt.subplots(
        1, n_cols, figsize=(2.05 * n_cols, 2.2), sharey=False,
    )
    if n_cols == 1: axes = [axes]

    for ax, (graph_label, data) in zip(axes, panels):
        for ratio, color in sorted(RATIO_COLORS.items()):
            rows = data.get((ALGO_KEY, ratio), [])
            if not rows: continue
            xs = [r[0] for r in rows]
            ys = [r[2] / 1024.0 if mem else r[1] for r in rows]
            ax.plot(xs, ys, color=color, marker="o", markersize=3.5,
                    linewidth=1.5, label=f"{int(ratio*100)}%")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.tick_params(axis="both", which="major", labelsize=8.5)
        ax.tick_params(axis="both", which="minor", labelsize=0)
        # Bold the x-axis tick labels per user request.
        for lbl in ax.get_xticklabels():
            lbl.set_fontweight("bold")
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.set_title(graph_label, fontsize=10, fontweight="bold")
        ax.set_xlabel(r"$\boldsymbol{s}$", fontsize=10, fontweight="bold")
    axes[0].set_ylabel(ylabel, fontsize=9)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(RATIO_COLORS),
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, 1.05),
               title="vertex sample ratio", title_fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot("time_ms",   "wall-clock time (ms)", OUT / "fig_scalability_time.pdf",
         mem=False)
    plot("memory_kB", "peak RSS (MB)",        OUT / "fig_scalability_mem.pdf",
         mem=True)
