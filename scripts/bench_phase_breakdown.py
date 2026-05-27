#!/usr/bin/env python3
"""Run degeneracy_cliques with PIVOTER_BREAKDOWN_LOG set; collect per-phase
time + RSS for Exp-4 (time breakdown) and Exp-5 (memory breakdown).

Output: one TSV containing
   meta  phase  duration_ms  rss_kb  delta_rss_kb  component_bytes
where `meta` encodes (graph, r, s, algo, run_id).

Reads the C++ PhaseLogger's TSV format produced by the existing
phaseStart()/phaseMark()/phaseDump() calls in main.cpp + RegND.

Usage:
    python3 scripts/bench_phase_breakdown.py \
        --graphs com-dblp web-it-2004 ca-HepPh ca-CondMat \
        --rs 3 \
        --ss 4 6 8 10 \
        --timeout 1800 \
        --out /data/wenqianz/phase_breakdown.tsv
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

BIN = "./build/bin/degeneracy_cliques"
PHASES_KEEP = ["loadAndSort", "buildSDCT", "preMutation", "prepareGraph",
               "MCEnum", "Index", "Support", "Peel", "dispatch_total"]


def run_one(graph: str, r: int, s: int, timeout_s: int,
            tsv_path: Path, run_id: int) -> dict:
    edges = f"graphs/{graph}.edges"
    if not os.path.exists(edges):
        return {"error": f"missing edges file: {edges}"}

    env = os.environ.copy()
    env["PIVOTER_RUN_REGION_V3LM"] = "1"
    env["PIVOTER_VSAFE_CLOUD"] = "1"
    env["PIVOTER_BREAKDOWN_LOG"] = str(tsv_path)
    env["PIVOTER_BREAKDOWN_META"] = f"{graph},{r},{s},V3LM,{run_id}"

    cmd = [BIN, edges, str(r), str(s), "degen"]
    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd, env=env, capture_output=True, text=True,
            timeout=timeout_s, check=False)
    except subprocess.TimeoutExpired:
        return {"error": f"timeout after {timeout_s}s"}
    dt = time.time() - t0
    rc = proc.returncode
    return {"wall_s": dt, "returncode": rc,
            "stdout_tail": proc.stdout[-400:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphs", nargs="+", required=True)
    ap.add_argument("--rs", type=int, nargs="+", required=True)
    ap.add_argument("--ss", type=int, nargs="+", required=True)
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--out", required=True,
                    help="path to TSV (PhaseLogger appends to it)")
    args = ap.parse_args()

    tsv_path = Path(args.out).resolve()
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    if tsv_path.exists():
        bak = tsv_path.with_suffix(tsv_path.suffix + f".bak.{int(time.time())}")
        tsv_path.rename(bak)
        print(f"[note] existing {tsv_path} renamed to {bak}", file=sys.stderr)

    run_id = 0
    n_done = 0; n_err = 0
    for g in args.graphs:
        for r in args.rs:
            for s in args.ss:
                if r >= s:
                    continue
                run_id += 1
                print(f"  [{g} r={r} s={s}] running... ", end="", flush=True)
                res = run_one(g, r, s, args.timeout, tsv_path, run_id)
                if res.get("error"):
                    print(f"ERR: {res['error']}")
                    n_err += 1
                else:
                    print(f"rc={res['returncode']} wall={res['wall_s']:.1f}s")
                    n_done += 1
    print(f"\ndone: {n_done} OK, {n_err} err, TSV at {tsv_path}")


if __name__ == "__main__":
    main()
