"""
End-to-end experiment figure in Nuclear-CD-TODS style.

Per-graph subplot: one subplot per graph, each plotting end-to-end
wall-clock time (or peak RSS) vs s, with two curves: SPIN-star (solid
blue, csv-key 'Pure') and CND (dashed red, csv-key 'REF_R1').  Output:
    fig_exp_time.pdf  -- time across all graphs (one row of subplots, may wrap)
    fig_exp_mem.pdf   -- peak RSS across all graphs (same layout)

Reads paper_data/01_main_benchmark_v3.csv (single source of truth for
r=1, three runs per cell; expects columns graph,r,s,algorithm,
wall_ms,...,time_max_rss_kB,...).  Filters status=OK; median across
runs per (graph, s, algorithm).

Type 42 fonts (VLDB/ACM compliant).
"""
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import csv
import math
from collections import defaultdict
from pathlib import Path

OUT  = Path(__file__).parent
CSV  = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data/01_main_benchmark_v3.csv")

# Graph display order + short labels (Nuclear-CD-style abbreviations).
GRAPH_ORDER = [
    ("com-dblp",                 "DBLP"),
    ("com-amazon.ungraph",       "AM"),
    ("twitter_combined",         "TW"),
    ("web-Stanford",             "STF"),
    ("web-Google",               "GO"),
    ("com-youtube",              "YT"),
    ("web-it-2004",              "WI"),
    ("wiki-Talk",                "WK"),
    ("soc-pokec-relationships",  "PO"),
    ("com-orkut",                "OR"),
]

ALGOS = [
    # csv-key,  paper-label,                  colour,   linestyle, marker, lw
    ("Pure",    r"SPIN$^{\star}$",            "#1f3a8a", "-",  "o", 1.4),
    ("REF_R1",  r"CND",                       "#c0392b", "--", "s", 1.2),
]


def load_data():
    """Returns {(graph, algo): [(s, time_ms, mem_kb), ...] sorted (median per s)}."""
    raw = defaultdict(lambda: defaultdict(list))
    if not CSV.exists():
        raise SystemExit(f"missing CSV: {CSV}")
    with CSV.open() as f:
        for r in csv.DictReader(f):
            if int(r.get("r", "0")) != 1: continue
            if r.get("status") != "OK":   continue
            try:
                s   = int(r["s"])
                # End-to-end time = wall_ms (the whole program: load, sort,
                # build, peel, output).  This is what the prose's "end-to-end"
                # claim refers to.
                t   = float(r["wall_ms"])
                m   = float(r["time_max_rss_kB"])
            except (ValueError, KeyError):
                continue
            raw[(r["graph"], r["algorithm"])][s].append((t, m))
    data = defaultdict(list)
    import statistics
    for key, per_s in raw.items():
        for s in sorted(per_s):
            ts = [x[0] for x in per_s[s]]
            ms = [x[1] for x in per_s[s]]
            data[key].append((s, statistics.median(ts), statistics.median(ms)))
    return data


def make_grid(metric, ylabel, out_pdf, log_y=True, mem=False):
    data = load_data()
    graphs = [(stem, lbl) for stem, lbl in GRAPH_ORDER
              if any((stem, a) in data for a, *_ in ALGOS)]
    n = len(graphs)
    # Lay out 5 per row.
    ncols = min(5, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 2.4, nrows * 1.2),
                             sharey=False)
    axes = axes.flatten() if nrows > 1 or ncols > 1 else [axes]

    for ax in axes[n:]:
        ax.set_visible(False)

    for idx, (stem, label) in enumerate(graphs):
        ax = axes[idx]
        for algo, alabel, color, ls, marker, lw in ALGOS:
            rows = data.get((stem, algo), [])
            if not rows: continue
            xs = [r[0] for r in rows]
            ys = [r[2] / 1024.0 if mem else r[1] for r in rows]
            ax.plot(xs, ys, color=color, linestyle=ls, marker=marker,
                    markersize=2.4, linewidth=lw, label=alabel,
                    markevery=max(1, len(xs)//12))
        ax.set_title(label, fontsize=8.5, pad=2)
        if log_y: ax.set_yscale("log")
        ax.set_xscale("log")
        ax.tick_params(axis="both", which="major", labelsize=6.5, pad=1)
        ax.tick_params(axis="both", which="minor", labelsize=0)
        ax.grid(True, which="major", alpha=0.25, linestyle=":")
        if idx % ncols == 0:
            ax.set_ylabel(ylabel, fontsize=7.5)
        if idx >= n - ncols:
            ax.set_xlabel(r"$s$", fontsize=7.5)

    # One shared legend at the top.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(ALGOS),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.96], h_pad=0.6, w_pad=0.4)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}  ({nrows}x{ncols} grid, {n} graphs)")


if __name__ == "__main__":
    make_grid("time_ms",   "wall-clock time (ms)", OUT / "fig_exp_time.pdf",
              log_y=True, mem=False)
    make_grid("memory_kB", "peak RSS (MB)",        OUT / "fig_exp_mem.pdf",
              log_y=True, mem=True)
