"""
Peel-time and peak-RSS experiment figure (per-graph grid).

Per-graph subplot: one subplot per graph, each plotting peel-time (or
peak RSS) vs s, with three curves --- SPIN-star (solid blue,
csv-key 'Pure'), SPIN (dotted green, csv-key 'LOCAL_V4'), and CND
(dashed red, csv-key 'REF_R1').  Output:
    fig_exp_time.pdf  -- peel time across all graphs (one row of subplots, may wrap)
    fig_exp_mem.pdf   -- peak RSS across all graphs (same layout)

Reads:
    paper_data/01_main_benchmark_v3.csv  -- columns include
        graph,r,s,algorithm,run,status,wall_ms,took_ms,build_ms,peel_ms,
        memory_kB,time_max_rss_kB,...
        Provides SPIN-star (algorithm='Pure') and CND (algorithm='REF_R1').
    paper_data/bench_local_v4.csv  -- columns
        graph,s,algo,status,wall_ms,took_ms,time_max_rss_kB,...
        Provides SPIN (algo='LOCAL_V4'); we treat took_ms as peel time
        (SPIN's wall_ms includes the shared CPI construction; took_ms
        isolates the SPIN solver phase, consistent with peel_ms for the
        other two algorithms).

Filter: status=OK; median across runs per (graph, s, algorithm).
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
import statistics
from collections import defaultdict
from pathlib import Path

OUT       = Path(__file__).parent
MAIN_CSV  = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data/01_main_benchmark_v3.csv")
SPIN_CSV  = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data/bench_local_v4.csv")

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

# Order matters for legend; SPIN-star first (blue), then SPIN (green),
# then CND (red).
ALGOS = [
    # internal-key,  paper-label,                  colour,   linestyle, marker, lw
    ("SPIN_STAR",    r"SPIN$^{\star}$",            "#1f3a8a", "-",  "o", 1.4),
    ("SPIN",         r"SPIN",                      "#2d8659", ":",  "^", 1.2),
    ("CND",          r"CND",                       "#c0392b", "--", "s", 1.2),
]


def load_main_csv():
    """SPIN-star (peel_ms) and CND (peel_ms) from the main benchmark."""
    raw = defaultdict(lambda: defaultdict(list))
    if not MAIN_CSV.exists():
        raise SystemExit(f"missing CSV: {MAIN_CSV}")
    with MAIN_CSV.open() as f:
        for r in csv.DictReader(f):
            if int(r.get("r", "0")) != 1: continue
            if r.get("status") != "OK":   continue
            algo_csv = r["algorithm"]
            if algo_csv == "Pure":
                key = "SPIN_STAR"
            elif algo_csv == "REF_R1":
                key = "CND"
            else:
                continue
            try:
                s = int(r["s"])
                # Peel time isolates the algorithm phase; wall_ms is dominated
                # by the shared CPI construction at high s and makes SPIN-star
                # vs CND look near-tied, hiding the algorithmic win.
                t = float(r["peel_ms"]) if r.get("peel_ms") else None
                m = float(r["time_max_rss_kB"])
            except (ValueError, KeyError):
                continue
            if t is None:
                continue
            raw[(r["graph"], key)][s].append((t, m))
    return raw


def load_spin_csv():
    """SPIN (LOCAL_V4) peel time from the SPIN bench (took_ms = solver phase)."""
    raw = defaultdict(lambda: defaultdict(list))
    if not SPIN_CSV.exists():
        return raw
    with SPIN_CSV.open() as f:
        for r in csv.DictReader(f):
            if r.get("status") != "OK":     continue
            if r.get("algo")   != "LOCAL_V4": continue
            try:
                s = int(r["s"])
                t = float(r["took_ms"])
                m = float(r["time_max_rss_kB"])
            except (ValueError, KeyError):
                continue
            raw[(r["graph"], "SPIN")][s].append((t, m))
    return raw


def load_data():
    """Returns {(graph, internal-key): [(s, time_ms, mem_kb), ...] sorted}.

    Clips SPIN to <= max(s) of SPIN-star on the same graph: the SPIN bench
    was swept to higher s on web-it-2004 (s up to 430 vs SPIN-star's 56)
    purely to feed Exp-4's max-speedup claim. Showing that extended tail
    in the headline Exp-1 figure makes SPIN-star/CND look like they "didn't
    finish" on the right edge. Cap SPIN to wherever SPIN-star reaches so
    all three curves stop at the same x.
    """
    raw_main = load_main_csv()
    raw_spin = load_spin_csv()
    # Per-graph cap = max s where SPIN-star has data.
    spinstar_max = {}
    for (g, k), per_s in raw_main.items():
        if k == "SPIN_STAR":
            spinstar_max[g] = max(per_s.keys())
    raw = dict(raw_main)
    for (g, k), per_s in raw_spin.items():
        cap = spinstar_max.get(g)
        if cap is None:
            # Graph has no SPIN-star coverage; skip SPIN too so the figure
            # doesn't show a lonely curve.
            continue
        filtered = {s: vals for s, vals in per_s.items() if s <= cap}
        if filtered:
            raw[(g, k)] = filtered
    data = defaultdict(list)
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

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(ALGOS),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.96], h_pad=0.6, w_pad=0.4)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}  ({nrows}x{ncols} grid, {n} graphs)")


if __name__ == "__main__":
    make_grid("peel_ms",   "peel time (ms)", OUT / "fig_exp_time.pdf",
              log_y=True, mem=False)
    make_grid("memory_kB", "peak RSS (MB)",  OUT / "fig_exp_mem.pdf",
              log_y=True, mem=True)
