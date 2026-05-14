"""
Stress-test figure: synthetic dense graphs (clique injection, N=1000),
density p ∈ {0.01..0.20}, time and RSS as functions of p.

Layout: 1 row × N subplots, one per s value.  Each subplot: time-vs-p
with two curves (Ours blue solid, SOTA red dashed).  A separate figure
shows memory-vs-p with the same layout.

Reads: paper_data/stress.csv with columns
    n,density,m,s,algorithm,time_ms,memory_kB,status
"""
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import csv
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).parent
CSV = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data/stress.csv")

ALGOS = [
    ("Ours_ST", "Ours", "#1f3a8a", "-",  "o", 1.4),
    ("REF_R1",  "SOTA", "#c0392b", "--", "s", 1.2),
]


def load():
    """Returns {(algo, s): [(p, time_ms, mem_kb), ...] sorted by p}."""
    rows = list(csv.DictReader(open(CSV)))
    data = defaultdict(list)
    for r in rows:
        if r.get("status") != "OK": continue
        try:
            p = float(r["density"])
            s = int(r["s"])
            t = float(r["time_ms"])
            m = float(r["memory_kB"])
        except (KeyError, ValueError):
            continue
        data[(r["algorithm"], s)].append((p, t, m))
    for k in data: data[k].sort()
    return data


def plot(metric, ylabel, out_pdf, mem=False):
    data = load()
    s_values = sorted({k[1] for k in data})
    n = len(s_values)
    fig, axes = plt.subplots(1, n, figsize=(n * 2.4, 1.9), sharey=False)
    if n == 1: axes = [axes]
    for ax, s in zip(axes, s_values):
        for algo, label, color, ls, marker, lw in ALGOS:
            rows = data.get((algo, s), [])
            if not rows: continue
            xs = [r[0] for r in rows]
            ys = [r[2] / 1024.0 if mem else r[1] for r in rows]
            ax.plot(xs, ys, color=color, linestyle=ls, marker=marker,
                    markersize=3.5, linewidth=lw, label=label)
        ax.set_title(f"$s{{=}}{s}$", fontsize=9)
        ax.set_xlabel("density $p$", fontsize=8)
        ax.set_yscale("log")
        ax.tick_params(axis="both", which="major", labelsize=7)
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
    axes[0].set_ylabel(ylabel, fontsize=8.5)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(ALGOS),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.05))
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot("time_ms", "wall-clock time (ms)", OUT / "fig_stress_time.pdf",
         mem=False)
    plot("memory_kB", "peak RSS (MB)",     OUT / "fig_stress_mem.pdf",
         mem=True)
