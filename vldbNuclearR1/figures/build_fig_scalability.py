"""
Scalability figure: time / RSS vs s, one curve per edge subsample ratio.

Two-subplot layout:
    Left:  Ours wall-clock vs s, 5 curves (20%, 40%, 60%, 80%, 100% edges)
    Right: SOTA wall-clock vs s, same 5 curves

The curves are expected to stack cleanly in proportion to the sample
ratio (near-linear scaling in edge count).

Reads: paper_data/scalability.csv with columns
    graph,base_edges,ratio,kept_edges,s,algorithm,run,time_ms,memory_kB,status
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

OUT = Path(__file__).parent
CSV = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data/scalability_com-dblp.csv")

ALGOS = [("Ours_ST", r"SPIN$^\star$ (ours)"), ("REF_R1", "CND")]

# 5 ratios -> 5 colours from a perceptually ordered palette.
RATIO_COLORS = {
    0.2: "#fee08b",
    0.4: "#fdae61",
    0.6: "#f46d43",
    0.8: "#d73027",
    1.0: "#a50026",
}


def load():
    """Returns {(algo, ratio): [(s, time_ms, mem_kb), ...] sorted}."""
    rows = list(csv.DictReader(open(CSV)))
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
    data = load()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True)
    for ax, (algo_key, algo_label) in zip(axes, ALGOS):
        for ratio, color in sorted(RATIO_COLORS.items()):
            rows = data.get((algo_key, ratio), [])
            if not rows: continue
            xs = [r[0] for r in rows]
            ys = [r[2] / 1024.0 if mem else r[1] for r in rows]
            ax.plot(xs, ys, color=color, marker="o", markersize=3.0,
                    linewidth=1.3, label=f"{int(ratio*100)}%")
        ax.set_title(algo_label, fontsize=10)
        ax.set_xlabel("$s$", fontsize=9)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.tick_params(axis="both", which="major", labelsize=7.5)
        ax.tick_params(axis="both", which="minor", labelsize=0)
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
    axes[0].set_ylabel(ylabel, fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(RATIO_COLORS),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.05),
               title="vertex sample ratio", title_fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot("time_ms", "wall-clock time (ms)", OUT / "fig_scalability_time.pdf",
         mem=False)
    plot("memory_kB", "peak RSS (MB)",     OUT / "fig_scalability_mem.pdf",
         mem=True)
