"""Parallel CPI construction figure: thread scaling per graph.

1 row x N cols.  Each panel shows wall-clock build time vs thread count
on a log-log axis; one line per s value tested for that graph.  Lines
extend to T=64 (the deepest sweep run); the slope on a log-log axis is
the speedup exponent.

Reads paper_data/bench_par_sdct.csv (status==OK rows).
"""
import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

OUT      = Path(__file__).parent
DATA_CSV = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data/bench_par_sdct.csv")

# Graphs selected for the figure (smallest-to-largest by edge count).
# dblp-core30 (overhead-dominated) and com-amazon are omitted: their
# build is too short for thread scaling to be informative.
GRAPHS = [
    ("com-dblp",                 "com-dblp"),
    ("twitter_combined",         "twitter"),
    ("com-youtube",              "com-youtube"),
    ("web-Google",               "web-Google"),
    ("web-Stanford",             "web-Stanford"),
    ("soc-pokec-relationships",  "soc-pokec"),
    ("web-it-2004",              "web-it-2004"),
]

# Pick at most two s values per graph: low and high, to keep panels readable.
S_PICKS = {
    "com-dblp":                 [3, 15],
    "twitter_combined":         [3, 12],
    "com-youtube":              [3, 16],
    "web-Google":               [3, 20],
    "web-Stanford":             [3, 60],
    "soc-pokec-relationships":  [3, 8],
    "web-it-2004":              [3, 400],
}

LOW_COLOR  = "#1f78b4"
HIGH_COLOR = "#e31a1c"


def load():
    """Returns {(graph, s, T): median build_ms}."""
    bucket = defaultdict(list)
    with open(DATA_CSV) as f:
        for r in csv.DictReader(f):
            if r.get("status") != "OK": continue
            try:
                T = int(r["threads"])
                build = float(r["build_ms"])
            except (ValueError, KeyError):
                continue
            if build <= 0: continue
            bucket[(r["graph"], int(r["s"]), T)].append(build)
    return {k: statistics.median(v) for k, v in bucket.items()}


def plot(out_pdf):
    data = load()
    n = len(GRAPHS)
    fig, axes = plt.subplots(1, n, figsize=(n * 2.0, 2.0), sharex=False)
    if n == 1: axes = [axes]

    for col, (stem, label) in enumerate(GRAPHS):
        ax = axes[col]
        for i, s in enumerate(S_PICKS[stem]):
            color = LOW_COLOR if i == 0 else HIGH_COLOR
            marker = "o" if i == 0 else "s"
            pts = sorted([(T, ms) for (g, ss, T), ms in data.items()
                          if g == stem and ss == s])
            if not pts: continue
            xs = [T for T, _ in pts]
            ys = [ms for _, ms in pts]
            ax.plot(xs, ys, color=color, marker=marker, markersize=4.0,
                    linewidth=1.4, label=f"$s{{=}}{s}$")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.tick_params(axis="both", which="major", labelsize=8.0)
        ax.tick_params(axis="both", which="minor", labelsize=0)
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.set_title(label, fontsize=9.5, fontweight="bold")
        ax.set_xlabel("threads", fontsize=9)
        ax.legend(fontsize=8, frameon=False, loc="lower left")
        ax.set_xticks([1, 4, 16, 64])
        if col == 0:
            ax.set_ylabel("build time (ms)", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot(OUT / "fig_par.pdf")
