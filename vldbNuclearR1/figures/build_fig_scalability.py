"""Combined scalability figure: SPIN★ time (row 1) and peak RSS (row 2) vs s
on five graphs, one curve per vertex-induced subsample ratio.

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

GRAPHS = [
    ("com-dblp",                 "com-dblp"),
    ("com-youtube",              "com-youtube"),
    ("web-Stanford",             "web-Stanford"),
    ("web-Google",               "web-Google"),
    ("soc-pokec-relationships",  "soc-pokec"),
]

ALGO_KEY = "Ours_ST"

RATIO_COLORS = {
    0.2: "#fee08b",
    0.4: "#fdae61",
    0.6: "#f46d43",
    0.8: "#d73027",
    1.0: "#a50026",
}


def load(csv_path):
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


def plot_combined(out_pdf):
    panels = []
    for stem, label in GRAPHS:
        d = load(DATA_DIR / f"scalability_{stem}.csv")
        if d:
            panels.append((label, d))
    if not panels: return

    n_cols = len(panels)
    fig, axes = plt.subplots(
        2, n_cols, figsize=(2.05 * n_cols, 4.0), sharey=False,
    )
    if n_cols == 1: axes = [[axes[0]], [axes[1]]]

    for row, (mem, ylabel) in enumerate([(False, "wall time (ms)"),
                                          (True,  "peak RSS (MB)")]):
        for col, (graph_label, data) in enumerate(panels):
            ax = axes[row][col]
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
            for lbl in ax.get_xticklabels():
                lbl.set_fontweight("bold")
            ax.grid(True, which="major", alpha=0.25, linestyle=":")
            for sp in ("top", "right"): ax.spines[sp].set_visible(False)
            if row == 0:
                ax.set_title(graph_label, fontsize=10, fontweight="bold")
            if row == 1:
                ax.set_xlabel(r"$\boldsymbol{s}$", fontsize=10, fontweight="bold")
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=9)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(RATIO_COLORS),
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, 1.03),
               title="vertex sample ratio", title_fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot_combined(OUT / "fig_scalability.pdf")
