"""Stress-test figure: synthetic dense graphs (clique injection, N=1000)
across density p and clique size s.

Emits a single fig_stress.pdf in a 2-row x n-col layout (time top,
memory bottom).  Each curve extends to the highest density at which the
algorithm still completes; cells with status TIMEOUT or FAIL are simply
dropped, so the line naturally ends at the limit of that algorithm.

Reads paper_data/stress.csv with columns
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
    ("Ours_ST", r"SPIN$^{\star}$",          "#1f3a8a", "-",  "o", 1.5),
    ("REF_R1",  r"CND",                     "#c0392b", "--", "s", 1.3),
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


def plot(out_pdf):
    data = load()
    s_values = sorted({k[1] for k in data})
    n = len(s_values)
    fig, axes = plt.subplots(2, n, figsize=(n * 2.4, 3.4),
                             sharex=False, sharey=False)
    if n == 1: axes = [[axes[0]], [axes[1]]]

    for row, (mem, ylabel) in enumerate([
        (False, "time (ms)"),
        (True,  "memory (MB)"),
    ]):
        for col, s in enumerate(s_values):
            ax = axes[row][col]
            for algo, label, color, ls, marker, lw in ALGOS:
                rows = data.get((algo, s), [])
                if not rows: continue
                xs = [r[0] for r in rows]
                ys = [r[2] / 1024.0 if mem else r[1] for r in rows]
                ax.plot(xs, ys, color=color, linestyle=ls, marker=marker,
                        markersize=4.0, linewidth=lw, label=label)
            ax.set_yscale("log")
            ax.tick_params(axis="both", which="major", labelsize=8.0)
            ax.tick_params(axis="both", which="minor", labelsize=0)
            ax.grid(True, which="major", alpha=0.25, linestyle=":")
            for sp in ("top", "right"): ax.spines[sp].set_visible(False)
            if row == 0:
                ax.set_title(f"$s{{=}}{s}$", fontsize=10, fontweight="bold")
            if row == 1:
                ax.set_xlabel("density $p$", fontsize=10)
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=9.5)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(ALGOS),
               fontsize=9.5, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot(OUT / "fig_stress.pdf")
