#!/usr/bin/env python3
"""
Time + memory breakdown experiment for R=1 (1,s)-nucleus decomposition.

For each (graph, s) configuration the script runs *both* algorithms
N times and aggregates the per-phase time / RSS / component-byte trace
emitted by src/PhaseLogger.h.  The two algorithms are:

    * ours    -- our tree-free centralized peel  (PIVOTER_RUN_ST_V2=1)
    * ref     -- the SOTA mutable-tree baseline   (no env flag)

The C++ binary writes one TSV row per phase to a path supplied via
PIVOTER_BREAKDOWN_LOG; we encode (graph, r, s, algo, run) in
PIVOTER_BREAKDOWN_META so rows can be aggregated later.

Outputs (under --out):
    breakdown_raw.tsv           -- every phase row from every run
    breakdown_median.csv        -- median per (graph, s, algo, phase)
    breakdown_summary.csv       -- one-line-per-(graph,s,algo) wall-clock + peak RSS
    breakdown_time_<graph>.pdf  -- stacked-bar time decomposition (Ours vs Ref) over s
    breakdown_mem_<graph>.pdf   -- stacked-bar memory decomposition (Ours vs Ref) over s

Usage:
    python3 scripts/run_breakdown.py \
        --bin   ./build/bin/degeneracy_cliques \
        --graph-dir /path/to/graphs \
        --graphs com-youtube web-Stanford \
        --s-list 5 10 15 20 25 \
        --runs 3 \
        --out  results/breakdown
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from collections import defaultdict
from statistics import median

# ---------------------------------------------------------------------------
# Raise stack limit so deep BK recursion (degeneracy ~400 on web-it-2004 at
# large s) doesn't segfault.  Modeled on bench_v3_all.py.
# ---------------------------------------------------------------------------
try:
    import resource
    _soft, _hard = resource.getrlimit(resource.RLIMIT_STACK)
    if _hard == resource.RLIM_INFINITY:
        _target = resource.RLIM_INFINITY
    else:
        _target = max(_soft, _hard - 4096)   # macOS rejects target == hard
    if _target > _soft:
        resource.setrlimit(resource.RLIMIT_STACK, (_target, _hard))
    _new_soft = resource.getrlimit(resource.RLIMIT_STACK)[0]
    def _fmt(x):
        return ("unlimited" if x == resource.RLIM_INFINITY
                else f"{x/1024/1024:.1f}MB")
    print(f"[stack] RLIMIT_STACK: {_fmt(_soft)} -> {_fmt(_new_soft)}",
          flush=True)
except Exception as e:
    print(f"[stack] WARNING: failed to raise stack limit: {e}", flush=True)

# Algorithm-key -> dict of (env vars to set) and label
ALGOS = {
    "ours": {
        "label":  "Ours (tree-free)",
        "env":    {"PIVOTER_RUN_ST_V2": "1"},
        # Phases that should be visualised as the algorithm's own time stack.
        # Order matters in the stacked plot.
        "phases_time": [
            "loadAndSort",        # graph load + degeneracy sort
            "STV2_SDCT_walk",     # Augmented SDCT walk + initial support fill
            "STV2_CSR_build",     # build dual CSR from COO
            "preMutation",        # other pre-mutation work (small)
            "prepareGraph",       # beSingleEdge etc.
            "STV2_peel_init",     # bucket array setup
            "STV2_peel_loop",     # peeling main loop
        ],
        "phases_mem":  [          # phases that contribute component_bytes for memory
            "STV2_SDCT_walk",
            "STV2_CSR_build",
            "STV2_peel_init",
        ],
    },
    "ref": {
        "label":  "Ref (SOTA mutable tree)",
        "env":    {},             # no PIVOTER_RUN_* flag → defaults dispatch to ref
        "phases_time": [
            "loadAndSort",
            "buildSDCT",          # SDCT tree build (the ref's Phase 2)
            "preMutation",
            "prepareGraph",
            "REF_initSupports",   # initial support computation
            "REF_heapBuild",      # heap + handle vector + leafRm structures
            "REF_peel_loop",      # peeling main loop
        ],
        "phases_mem":  [
            "buildSDCT",
            "REF_initSupports",
            "REF_heapBuild",
        ],
    },
}


# ---------------------------------------------------------------------------
# Run subprocess for a single (graph, s, algo, run) and append to raw_tsv.
# ---------------------------------------------------------------------------

def run_one(bin_path: Path, graph_path: Path, s: int, algo_key: str,
            run_id: int, raw_tsv: Path, timeout: float) -> tuple[bool, float, float]:
    """Run the binary once.  Returns (ok, wall_ms, peak_rss_kb)."""
    algo = ALGOS[algo_key]
    env = os.environ.copy()
    env.update(algo["env"])
    env["PIVOTER_BREAKDOWN_LOG"] = str(raw_tsv)
    env["PIVOTER_BREAKDOWN_META"] = (
        f"graph={graph_path.stem},r=1,s={s},algo={algo_key},run={run_id}"
    )

    t0 = time.time()
    try:
        proc = subprocess.run(
            [str(bin_path), str(graph_path), "1", str(s)],
            env=env, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT after {timeout:.0f}s")
        return False, -1.0, -1.0

    wall_ms = (time.time() - t0) * 1000.0
    if proc.returncode != 0:
        print(f"  FAIL rc={proc.returncode}")
        print((proc.stderr or "")[-400:])
        return False, wall_ms, -1.0

    # Parse "Final Memory: NNNN kB" from stdout for an end-of-run RSS reading.
    peak_rss_kb = -1.0
    for line in proc.stdout.splitlines():
        if "Final Memory:" in line:
            try:
                peak_rss_kb = float(line.split(":")[-1].strip().split()[0])
            except Exception:
                pass
    return True, wall_ms, peak_rss_kb


# ---------------------------------------------------------------------------
# Aggregation -- read raw TSV, group by (graph, s, algo, phase), median across runs
# ---------------------------------------------------------------------------

def parse_meta(meta: str) -> dict:
    out = {}
    for kv in meta.split(","):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def aggregate(raw_tsv: Path) -> dict:
    """Returns {(graph, s, algo): {phase: {'time_ms': median, 'bytes': median, 'rss_kb': median}}}"""
    grouped = defaultdict(lambda: defaultdict(list))
    with raw_tsv.open() as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            meta_str, phase, dur, rss, drss, bytes_ = parts[:6]
            m = parse_meta(meta_str)
            key = (m.get("graph"), int(m.get("s", -1)), m.get("algo"))
            grouped[key][phase].append({
                "time_ms":  float(dur),
                "rss_kb":   float(rss),
                "delta_kb": float(drss),
                "bytes":    int(bytes_),
            })

    result = {}
    for key, phases in grouped.items():
        agg = {}
        for ph, runs in phases.items():
            agg[ph] = {
                "time_ms": median(r["time_ms"]  for r in runs),
                "rss_kb":  median(r["rss_kb"]   for r in runs),
                "delta_kb": median(r["delta_kb"] for r in runs),
                "bytes":   int(median(r["bytes"]   for r in runs)),
                "n_runs":  len(runs),
            }
        result[key] = agg
    return result


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_median(agg, path: Path):
    cols = ["graph", "s", "algo", "phase", "n_runs",
            "time_ms", "rss_kb", "delta_rss_kb", "component_bytes"]
    with path.open("w") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for (g, s, a), phases in sorted(agg.items()):
            for ph, v in phases.items():
                w.writerow([g, s, a, ph, v["n_runs"],
                            f"{v['time_ms']:.3f}",
                            f"{v['rss_kb']:.0f}",
                            f"{v['delta_kb']:.0f}",
                            v["bytes"]])
    print(f"wrote {path}")


def write_summary(agg, path: Path):
    """One row per (graph, s, algo) with total time + peak RSS."""
    cols = ["graph", "s", "algo", "total_time_ms", "peak_rss_kb",
            "load_ms", "build_ms", "peel_ms"]
    with path.open("w") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for (g, s, a), phases in sorted(agg.items()):
            total = sum(v["time_ms"] for v in phases.values())
            peak  = max(v["rss_kb"]  for v in phases.values()) if phases else 0
            load  = phases.get("loadAndSort", {}).get("time_ms", 0)
            if a == "ours":
                build = (phases.get("STV2_SDCT_walk", {}).get("time_ms", 0)
                         + phases.get("STV2_CSR_build", {}).get("time_ms", 0))
                peel  = (phases.get("STV2_peel_init", {}).get("time_ms", 0)
                         + phases.get("STV2_peel_loop", {}).get("time_ms", 0))
            else:
                build = (phases.get("buildSDCT", {}).get("time_ms", 0)
                         + phases.get("REF_initSupports", {}).get("time_ms", 0)
                         + phases.get("REF_heapBuild", {}).get("time_ms", 0))
                peel  = phases.get("REF_peel_loop", {}).get("time_ms", 0)
            w.writerow([g, s, a,
                        f"{total:.2f}", f"{peak:.0f}",
                        f"{load:.2f}", f"{build:.2f}", f"{peel:.2f}"])
    print(f"wrote {path}")


# ---------------------------------------------------------------------------
# Plots: stacked-bar comparing Ours vs Ref over s for each graph
# ---------------------------------------------------------------------------

def make_plots(agg, out_dir: Path, graphs: list[str], s_list: list[int]):
    try:
        import matplotlib
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"]  = 42
        matplotlib.rcParams["text.usetex"]  = False
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[plot] matplotlib not available; skipping figures")
        return

    palette = [
        "#3a5fcd",  # blue
        "#1f8a3a",  # green
        "#e8a44c",  # orange
        "#7a2bb0",  # purple
        "#b05020",  # rust
        "#888888",  # grey
        "#ce375a",  # pink
    ]

    for g in graphs:
        # ===== Time stacked plot =====
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True)
        for ax, algo_key in zip(axes, ["ours", "ref"]):
            algo = ALGOS[algo_key]
            x = np.arange(len(s_list))
            bottom = np.zeros(len(s_list))
            for i, ph in enumerate(algo["phases_time"]):
                heights = []
                for s in s_list:
                    v = agg.get((g, s, algo_key), {}).get(ph, {}).get("time_ms", 0)
                    heights.append(v)
                ax.bar(x, heights, bottom=bottom, label=ph,
                       color=palette[i % len(palette)],
                       edgecolor="black", linewidth=0.4, width=0.7)
                bottom += np.array(heights)
            ax.set_xticks(x)
            ax.set_xticklabels([str(s) for s in s_list])
            ax.set_xlabel("$s$")
            ax.set_title(algo["label"])
            ax.grid(axis="y", alpha=0.2, linestyle=":")
            ax.legend(fontsize=7, loc="upper left", frameon=False)
        axes[0].set_ylabel("wall-clock time (ms)")
        fig.suptitle(f"{g}: time breakdown by phase", fontsize=11)
        fig.tight_layout()
        path = out_dir / f"breakdown_time_{g}.pdf"
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {path}")

        # ===== Memory stacked plot =====
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True)
        for ax, algo_key in zip(axes, ["ours", "ref"]):
            algo = ALGOS[algo_key]
            x = np.arange(len(s_list))
            bottom = np.zeros(len(s_list))
            for i, ph in enumerate(algo["phases_mem"]):
                heights = []
                for s in s_list:
                    v = agg.get((g, s, algo_key), {}).get(ph, {}).get("bytes", 0)
                    heights.append(v / (1024 * 1024))   # MB
                ax.bar(x, heights, bottom=bottom, label=ph,
                       color=palette[i % len(palette)],
                       edgecolor="black", linewidth=0.4, width=0.7)
                bottom += np.array(heights)
            ax.set_xticks(x)
            ax.set_xticklabels([str(s) for s in s_list])
            ax.set_xlabel("$s$")
            ax.set_title(algo["label"])
            ax.grid(axis="y", alpha=0.2, linestyle=":")
            ax.legend(fontsize=7, loc="upper left", frameon=False)
        axes[0].set_ylabel("dominant data structures (MB)")
        fig.suptitle(f"{g}: memory breakdown by component", fontsize=11)
        fig.tight_layout()
        path = out_dir / f"breakdown_mem_{g}.pdf"
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin",       default="./build/bin/degeneracy_cliques",
                    help="Path to the degeneracy_cliques binary")
    ap.add_argument("--graph-dir", default=".",
                    help="Directory holding the .edges files")
    ap.add_argument("--graphs",    nargs="+", required=True,
                    help="Graph stem names; each will be looked up as <graph-dir>/<name>.edges")
    ap.add_argument("--s-list",    nargs="+", type=int, required=True,
                    help="Values of s to test")
    ap.add_argument("--algos",     nargs="+", default=["ours", "ref"],
                    choices=list(ALGOS.keys()),
                    help="Which algorithms to run (default: both)")
    ap.add_argument("--runs",      type=int, default=3,
                    help="Number of repetitions per (graph,s,algo) for median")
    ap.add_argument("--timeout",   type=float, default=1800.0,
                    help="Per-run timeout in seconds (default 1800)")
    ap.add_argument("--out",       default="results/breakdown",
                    help="Output directory")
    ap.add_argument("--skip-run",  action="store_true",
                    help="Skip subprocess execution; just re-aggregate from existing raw TSV")
    ap.add_argument("--resume",    action="store_true",
                    help="Resume by skipping (graph,s,algo,run) tuples already in the raw TSV")
    ap.add_argument("--fresh",     action="store_true",
                    help="Start fresh: truncate the raw TSV before running (default unless --resume)")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_tsv = out_dir / "breakdown_raw.tsv"

    bin_path = Path(args.bin)
    if not bin_path.exists():
        sys.exit(f"binary not found: {bin_path}")
    graph_dir = Path(args.graph_dir)

    # Resolve graph file extensions: prefer .edges, fall back to .edges.gz / .txt / .ungraph.txt
    def resolve(stem: str) -> Path | None:
        for ext in (".edges", ".edges.gz", ".ungraph.txt", ".txt"):
            p = graph_dir / (stem + ext)
            if p.exists(): return p
        return None

    # Read which (graph, s, algo, run) tuples are already done in the existing TSV.
    done_tuples = set()
    if args.resume and raw_tsv.exists():
        with raw_tsv.open() as f:
            f.readline()
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if not parts: continue
                m = parse_meta(parts[0])
                key = (m.get("graph"), m.get("s"), m.get("algo"), m.get("run"))
                done_tuples.add(key)
        print(f"[resume] {len(done_tuples)} (graph,s,algo,run) phase rows already in TSV; "
              f"will skip those exact tuples", flush=True)

    if not args.skip_run:
        # Truncate raw log unless --resume requested.
        if not args.resume and raw_tsv.exists():
            raw_tsv.unlink()

        # === run loop with skip-on-timeout-propagation per (graph, algo) ===
        n_total = len(args.graphs) * len(args.s_list) * len(args.algos) * args.runs
        n_done  = 0
        # Track for each (graph, algo) pair the smallest s that timed out;
        # everything s'>= that threshold gets skipped (they only get harder).
        timeout_floor = {}   # (graph, algo) -> int s
        for g in args.graphs:
            gp = resolve(g)
            if gp is None:
                print(f"[skip] graph {g}: no file in {graph_dir}", flush=True)
                continue
            for s in args.s_list:
                for algo_key in args.algos:
                    floor = timeout_floor.get((g, algo_key))
                    if floor is not None and s >= floor:
                        for run_id in range(args.runs):
                            n_done += 1
                            print(f"[{n_done}/{n_total}] {g} s={s} algo={algo_key}: "
                                  f"SKIP (earlier s={floor} timed out)", flush=True)
                        continue
                    for run_id in range(args.runs):
                        n_done += 1
                        tup = (g, str(s), algo_key, str(run_id))
                        if tup in done_tuples:
                            print(f"[{n_done}/{n_total}] {g} s={s} algo={algo_key} run={run_id}: "
                                  f"SKIP (already in TSV)", flush=True)
                            continue
                        print(f"[{n_done}/{n_total}] {g} s={s} algo={algo_key} run={run_id} ... ",
                              end="", flush=True)
                        ok, wall_ms, peak_rss = run_one(
                            bin_path, gp, s, algo_key, run_id, raw_tsv, args.timeout)
                        if ok:
                            print(f"OK  wall={wall_ms:.0f}ms  rss={peak_rss:.0f}kB", flush=True)
                        else:
                            # If this run timed out, propagate skip to larger s.
                            if wall_ms >= args.timeout * 1000 - 100:
                                timeout_floor[(g, algo_key)] = s
                                print(f"[skip-floor] {g} {algo_key} s>={s} will be skipped",
                                      flush=True)

    if not raw_tsv.exists():
        sys.exit(f"no raw TSV at {raw_tsv}")

    # === aggregate ===
    agg = aggregate(raw_tsv)
    write_median(agg,  out_dir / "breakdown_median.csv")
    write_summary(agg, out_dir / "breakdown_summary.csv")

    # === plots ===
    make_plots(agg, out_dir, args.graphs, args.s_list)


if __name__ == "__main__":
    main()
