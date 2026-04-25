#!/usr/bin/env python3
"""
Convert breakdown_median.csv from run_breakdown.py into LaTeX tables ready
to paste into §7 of the paper.

Emits two tables to stdout (or to --out):
    1. Time decomposition: per (graph, s) row, columns load / build / peel /
       total for both algorithms, plus ratio columns.
    2. Memory decomposition: per (graph, s) row, columns peak_RSS for both
       algorithms plus ratio.

Also emits a per-component memory table for one chosen graph that breaks
out the data-structure-level component_bytes (CSR vs hash, bucket vs heap,
tree vs none).
"""
import argparse, csv, sys
from collections import defaultdict
from pathlib import Path


def kb_to_mb(kb):
    try:
        return float(kb) / 1024.0
    except (ValueError, TypeError):
        return 0.0


def load_summary(path: Path):
    rows = list(csv.DictReader(open(path)))
    out = {}
    for r in rows:
        key = (r["graph"], int(r["s"]), r["algo"])
        out[key] = r
    return out


def load_median(path: Path):
    """Returns {(graph, s, algo): {phase: row}}."""
    rows = list(csv.DictReader(open(path)))
    out = defaultdict(dict)
    for r in rows:
        out[(r["graph"], int(r["s"]), r["algo"])][r["phase"]] = r
    return out


def emit_time_table(summary, graphs, out=sys.stdout):
    print(r"\begin{table}[t]", file=out)
    print(r"\centering", file=out)
    print(r"\caption{Per-phase wall-clock time (ms), median over 3 runs.}", file=out)
    print(r"\label{tab:bd-time}", file=out)
    print(r"\small", file=out)
    print(r"\begin{tabular}{llrrrrrr}", file=out)
    print(r"\toprule", file=out)
    print(r"\multirow{2}{*}{Graph} & \multirow{2}{*}{$s$} "
          r"& \multicolumn{3}{c}{Ours} & \multicolumn{3}{c}{Ref (SOTA)} \\", file=out)
    print(r"\cmidrule(l){3-5}\cmidrule(l){6-8}", file=out)
    print(r"& & load & build & peel & load & build & peel \\", file=out)
    print(r"\midrule", file=out)
    for g in graphs:
        ss = sorted({k[1] for k in summary if k[0] == g})
        for s in ss:
            ours = summary.get((g, s, "ours"), {})
            ref  = summary.get((g, s, "ref"),  {})
            row = (
                fr"\texttt{{{g}}} & {s} & "
                f"{float(ours.get('load_ms', 0)):.1f} & "
                f"{float(ours.get('build_ms', 0)):.1f} & "
                f"{float(ours.get('peel_ms', 0)):.2f} & "
                f"{float(ref.get('load_ms', 0)):.1f} & "
                f"{float(ref.get('build_ms', 0)):.1f} & "
                f"{float(ref.get('peel_ms', 0)):.2f} \\\\"
            )
            print(row, file=out)
        if g != graphs[-1]: print(r"\midrule", file=out)
    print(r"\bottomrule", file=out)
    print(r"\end{tabular}", file=out)
    print(r"\end{table}", file=out)


def emit_memory_table(summary, graphs, out=sys.stdout):
    print(file=out)
    print(r"\begin{table}[t]", file=out)
    print(r"\centering", file=out)
    print(r"\caption{Peak resident-set memory (MB), median over 3 runs.}", file=out)
    print(r"\label{tab:bd-mem}", file=out)
    print(r"\small", file=out)
    print(r"\begin{tabular}{llrrr}", file=out)
    print(r"\toprule", file=out)
    print(r"Graph & $s$ & Ours (MB) & Ref (MB) & Ratio \\", file=out)
    print(r"\midrule", file=out)
    for g in graphs:
        ss = sorted({k[1] for k in summary if k[0] == g})
        for s in ss:
            ours = summary.get((g, s, "ours"), {})
            ref  = summary.get((g, s, "ref"),  {})
            o_mb = kb_to_mb(ours.get("peak_rss_kb", 0))
            r_mb = kb_to_mb(ref.get("peak_rss_kb", 0))
            ratio = (r_mb / o_mb) if o_mb else 0
            print(fr"\texttt{{{g}}} & {s} & {o_mb:.0f} & {r_mb:.0f} & "
                  fr"{ratio:.2f}$\times$ \\", file=out)
        if g != graphs[-1]: print(r"\midrule", file=out)
    print(r"\bottomrule", file=out)
    print(r"\end{tabular}", file=out)
    print(r"\end{table}", file=out)


def emit_component_table(median, graph, out=sys.stdout):
    print(file=out)
    print(r"\begin{table}[t]", file=out)
    print(r"\centering", file=out)
    print(fr"\caption{{Per-component memory (MB) on \texttt{{{graph}}}, "
          r"median over 3 runs. Rows tagged with `--' do not exist for "
          r"the column's algorithm.}}", file=out)
    print(r"\label{tab:bd-component}", file=out)
    print(r"\small", file=out)
    keys_sorted = sorted({k for k in median if k[0] == graph}, key=lambda k: k[1])
    s_values = sorted({k[1] for k in keys_sorted})
    print(r"\begin{tabular}{l" + "r" * len(s_values) + r"}", file=out)
    print(r"\toprule", file=out)
    print(r"Component (algorithm) & " + " & ".join([f"$s{{=}}{s}$" for s in s_values])
          + r" \\", file=out)
    print(r"\midrule", file=out)

    # Map each (component, algo_key) to the phase name in the median CSV.
    rows = [
        ("CPI tree (Ref)",        "ref",  "buildSDCT"),
        ("Initial supports (Ref)", "ref", "REF_initSupports"),
        ("Heap + handles (Ref)", "ref",  "REF_heapBuild"),
        (r"\midrule", None, None),
        ("SDCT walk + COO (Ours)", "ours", "STV2_SDCT_walk"),
        ("Dual CSR + supports (Ours)", "ours", "STV2_CSR_build"),
        ("Bucket array (Ours)",  "ours", "STV2_peel_init"),
    ]
    for label, algo, phase in rows:
        if label == r"\midrule":
            print(r"\midrule", file=out)
            continue
        cells = []
        for s in s_values:
            row = median.get((graph, s, algo), {}).get(phase, {})
            mb = float(row.get("component_bytes", 0)) / (1024 * 1024)
            cells.append(f"{mb:.1f}" if mb > 0.05 else "--")
        print(label + " & " + " & ".join(cells) + r" \\", file=out)
    print(r"\bottomrule", file=out)
    print(r"\end{tabular}", file=out)
    print(r"\end{table}", file=out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", default="results/breakdown/breakdown_summary.csv")
    ap.add_argument("--median",  default="results/breakdown/breakdown_median.csv")
    ap.add_argument("--graphs",  nargs="+", default=None,
                    help="Graphs to include in the summary tables; default: all")
    ap.add_argument("--component-graph", default="com-youtube",
                    help="Graph to use for the per-component memory table")
    ap.add_argument("--out",     default=None, help="Output file (default stdout)")
    args = ap.parse_args()

    summary = load_summary(Path(args.summary))
    median  = load_median(Path(args.median))
    if args.graphs:
        graphs = args.graphs
    else:
        graphs = sorted({k[0] for k in summary})

    out = open(args.out, "w") if args.out else sys.stdout
    emit_time_table(summary, graphs, out)
    emit_memory_table(summary, graphs, out)
    emit_component_table(median, args.component_graph, out)
    if args.out: out.close()


if __name__ == "__main__":
    main()
