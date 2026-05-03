#!/usr/bin/env python3
"""
Summarise paper_data/bench_par_sdct.csv into a publishable speedup table.

Usage:
    python3 analyze_par_sdct.py              # auto-locates the CSV
    python3 analyze_par_sdct.py path.csv     # explicit path

Outputs:
- A per-(graph, s) table of speedup = median(build_ms_T1) / median(build_ms_T)
- Geometric-mean speedup across all (graph, s) cells, per T
- Best-T pick per (graph, s) — i.e. where speedup peaks before NUMA hits
- Honest report on cells with insufficient runs
"""
from __future__ import annotations
import csv, math, statistics, sys
from collections import defaultdict
from pathlib import Path


def load(csv_path: Path):
    rows = list(csv.DictReader(csv_path.open()))
    by_cell: dict = defaultdict(list)
    for r in rows:
        if r.get("status") != "OK" or not r.get("build_ms"): continue
        try:
            k = (r["graph"], int(r["s"]), int(r["threads"]))
            by_cell[k].append(int(r["build_ms"]))
        except (KeyError, ValueError):
            continue
    return by_cell


def gmean(xs):
    xs = [x for x in xs if x and x > 0]
    if not xs: return float("nan")
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def main():
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 \
               else Path("paper_data/bench_par_sdct.csv")
    if not csv_path.exists():
        print(f"ERROR: {csv_path} not found", file=sys.stderr); sys.exit(1)

    cells = load(csv_path)
    median = {k: int(statistics.median(v)) for k, v in cells.items()}
    counts = {k: len(v) for k, v in cells.items()}

    THREADS = [1, 2, 4, 8, 16, 24, 32, 48, 64]
    graphs_s = sorted({(g, s) for (g, s, _) in median.keys()})

    print("=" * 100)
    print(f"  bench_par_sdct.csv summary  ({csv_path})")
    print("=" * 100)
    print(f"  total cells: {len(median)}   (graph, s) pairs: {len(graphs_s)}")
    insufficient = [k for k, n in counts.items() if n < 3]
    if insufficient:
        print(f"  WARNING: {len(insufficient)} cells have <3 runs (still in progress?)")

    # Per-(graph, s) speedup table.
    print()
    print(f"{'graph':25s} {'s':>3s}   T=1(ms)   " + "  ".join(f"T={t:>2d}" for t in THREADS[1:]))
    print("-" * 100)
    by_T_speedups = defaultdict(list)
    best_T = []
    for g, s in graphs_s:
        t1 = median.get((g, s, 1))
        if not t1: continue
        line = f"{g:25s} {s:>3d}   {t1:>7d}   "
        sp_row = []
        for T in THREADS[1:]:
            tT = median.get((g, s, T))
            if tT:
                sp = t1 / tT
                line += f"{sp:>4.1f}x  "
                sp_row.append((T, sp))
                by_T_speedups[T].append(sp)
            else:
                line += "   --   "
                sp_row.append((T, None))
        # find peak T
        sp_only = [(T, sp) for T, sp in sp_row if sp is not None]
        if sp_only:
            bT, bS = max(sp_only, key=lambda x: x[1])
            best_T.append((g, s, bT, bS))
            line += f"  best=T{bT}({bS:.1f}x)"
        print(line)

    # Geometric mean across cells per T.
    print()
    print("=== Geometric-mean speedup across all (graph, s) cells, per T ===")
    for T in THREADS[1:]:
        vals = by_T_speedups[T]
        if vals:
            print(f"  T={T:>2d}: gmean={gmean(vals):>5.2f}x  ({len(vals)} cells, max={max(vals):.1f}x)")

    # Where does the peak land?
    print()
    print("=== Distribution of best-scaling T (peak before NUMA hits) ===")
    bestT_count = defaultdict(int)
    for _, _, T, _ in best_T:
        bestT_count[T] += 1
    for T in THREADS[1:]:
        if bestT_count[T]:
            print(f"  T={T:>2d}: {bestT_count[T]:>3d} cells peak here")

    # Top 10 highest absolute speedups.
    print()
    print("=== Top 10 highest speedups (graph, s, T) ===")
    triples = [(median[(g,s,1)] / median[(g,s,T)], g, s, T)
               for (g, s, T) in median if T != 1 and (g,s,1) in median]
    for sp, g, s, T in sorted(triples, reverse=True)[:10]:
        print(f"  {g:25s} s={s:>3d} T={T:>2d}: {sp:>5.1f}x  "
              f"({median[(g,s,1)]/1000:.1f}s -> {median[(g,s,T)]/1000:.2f}s)")

    print()
    print("=" * 100)


if __name__ == "__main__":
    main()
