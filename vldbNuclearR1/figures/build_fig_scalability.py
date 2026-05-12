"""
Scalability figure: time / RSS vs s, one curve per vertex-induced subsample
ratio, on multiple graphs.

Layout: rows = graphs, two columns (SPIN★ / CND). Single figure shows time
or RSS depending on `metric`. Reads
paper_data/scalability_{graph}.csv for each graph in GRAPHS.
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
    ("com-orkut",                "com-orkut"),
]

ALGOS = [("Ours_ST", r"SPIN$^\star$ (ours)"), ("REF_R1", "CND")]

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
    # Load all graphs first; drop any that has no CSV yet.
    panels = []
    for stem, label in GRAPHS:
        d = load(DATA_DIR / f"scalability_{stem}.csv")
        if d:
            panels.append((label, d))
    if not panels: return

    n_rows = len(panels)
    fig, axes = plt.subplots(
        n_rows, 2, figsize=(7.0, 1.55 * n_rows + 0.4),
        sharex=True, sharey="row",
    )
    if n_rows == 1:
        axes = [axes]  # normalise to 2D
    for row_idx, (graph_label, data) in enumerate(panels):
        for col_idx, (algo_key, algo_label) in enumerate(ALGOS):
            ax = axes[row_idx][col_idx]
            for ratio, color in sorted(RATIO_COLORS.items()):
                rows = data.get((algo_key, ratio), [])
                if not rows: continue
                xs = [r[0] for r in rows]
                ys = [r[2] / 1024.0 if mem else r[1] for r in rows]
                ax.plot(xs, ys, color=color, marker="o", markersize=3.0,
                        linewidth=1.3, label=f"{int(ratio*100)}%")
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
            ax.xaxis.set_minor_formatter(mticker.NullFormatter())
            ax.tick_params(axis="both", which="major", labelsize=7.5)
            ax.tick_params(axis="both", which="minor", labelsize=0)
            ax.grid(True, which="major", alpha=0.25, linestyle=":")
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            # Column titles only on the first row.
            if row_idx == 0:
                ax.set_title(algo_label, fontsize=10, fontweight="bold")
            # x labels only on bottom row.
            if row_idx == n_rows - 1:
                ax.set_xlabel(r"$s$", fontsize=9)
        # y label only on the left column; combine the graph label into it.
        axes[row_idx][0].set_ylabel(f"{graph_label}\n{ylabel}", fontsize=9)

    # Single shared legend above the figure.
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(RATIO_COLORS),
               fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, 1.0 + 0.5/(1.55*n_rows)),
               title="vertex sample ratio", title_fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 1 - 0.5/(1.55*n_rows)])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot("time_ms",   "wall-clock time (ms)", OUT / "fig_scalability_time.pdf",
         mem=False)
    plot("memory_kB", "peak RSS (MB)",        OUT / "fig_scalability_mem.pdf",
         mem=True)
