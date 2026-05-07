#!/usr/bin/env python3
"""Minimal serial bench identical in approach to test_mimic.py which proved
it survives. Drops the threading.Lock + function-def closures that the
fancier main() had — those (or something near them) have been correlating
with mystery SIGKILLs on tods2 during the past several hours of debugging.

Iterates the same job list as bench_r1_v3.main() and writes the same CSV.
Uses bench_r1_v3.run_one() directly so log paths and parse stay unified."""
from __future__ import annotations
import csv as _csv, os, sys, time
from pathlib import Path
sys.path.insert(0, '/home/wenqianz/nclique')
import bench_r1_v3 as b

OUTCSV = b.OUTCSV
LOGDIR = b.LOGDIR
GRAPHS = b.GRAPHS
ALGOS  = b.ALGOS
RUNS_PER_CFG = b.RUNS_PER_CFG
FIELDNAMES   = b.FIELDNAMES


def main():
    print(f"[mini] start {time.strftime('%F %T')}", flush=True)
    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    avail = []
    for g, smax in GRAPHS:
        if Path(f"./graphs/{g}.edges").exists():
            avail.append((g, smax))
    print(f"[mini] graphs: {[g for g, _ in avail]}", flush=True)

    counts = b.load_existing_counts(OUTCSV)
    n_done = sum(1 for n in counts.values() if n >= RUNS_PER_CFG)
    print(f"[mini] {n_done} (graph,s,algo) cells already complete", flush=True)

    jobs = []
    for g, smax in avail:
        for s in range(2, smax + 1):
            for algo, env_extra in ALGOS:
                already = counts.get((g, s, algo), 0)
                for run_idx in range(already, RUNS_PER_CFG):
                    jobs.append((g, s, algo, env_extra, run_idx))
    print(f"[mini] {len(jobs)} jobs to run", flush=True)

    skip_floor = {}
    for i, job in enumerate(jobs, 1):
        graph, s, algo, env_extra, run_idx = job
        if (graph, algo) in skip_floor and s >= skip_floor[(graph, algo)]:
            row = {"graph": graph, "r": 1, "s": s, "algorithm": algo,
                   "run": run_idx, "status": "SKIP_FLOOR"}
            b.append_row(OUTCSV, row)
            print(f"  [{i}/{len(jobs)}] {graph} s={s} {algo} r={run_idx} SKIP_FLOOR",
                  flush=True)
            continue
        result = b.run_one(graph, s, algo, env_extra, run_idx)
        row = {
            "graph": graph, "r": 1, "s": s, "algorithm": algo,
            "run": run_idx, "status": result["status"],
            "wall_ms": f"{result['wall_ms']:.1f}",
            "took_ms":  result.get("took_ms",  ""),
            "build_ms": result.get("build_ms", ""),
            "peel_ms":  result.get("peel_ms",  ""),
            "memory_kB": result.get("memory_kB", ""),
            **{k: result.get(k, "") for k in (
                "time_max_rss_kB","time_user_sec","time_sys_sec",
                "time_elapsed","time_pagefaults_major",
                "time_pagefaults_minor","time_voluntary_ctxt",
                "time_involuntary_ctxt","time_exit_status")},
        }
        b.append_row(OUTCSV, row)
        print(f"  [{i}/{len(jobs)}] {graph} s={s} {algo} r={run_idx} "
              f"{result['status']} wall={result['wall_ms']/1000:.1f}s",
              flush=True)
        if result["status"] == "TIMEOUT" or result["status"].startswith("FAIL(-9"):
            skip_floor[(graph, algo)] = s
            print(f"    -> skip-floor: {graph}/{algo} s>={s}", flush=True)
    print(f"[mini] done", flush=True)


if __name__ == "__main__":
    main()
