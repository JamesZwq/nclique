#!/usr/bin/env python3
"""
Experiment 1 — Local-iterative (H-index) study (Experimental.tex §6.3).

Goals:
  (A) T(s): iteration count of the parallel H-index algorithm vs s, on the
            same graphs as the main bench, single-thread.
  (B) Parallel-scaling: wall-clock of LOCAL_V4 at OMP_NUM_THREADS in
            {1,2,4,8,16,24} on com-dblp, com-youtube, web-Stanford.
            Plus ST_V3 (centralised peel) at the same configs as a baseline.

Phase A is parallel-outer (many jobs at once, each OMP=1).
Phase B is sequential-outer (one job at a time, each binding all cores).

Usage:
    nohup python3 bench_local_iterative.py tods2 > bench_local_iterative.log 2>&1 &
"""
from __future__ import annotations
import os, subprocess, sys, time
from pathlib import Path

from bench_lib import (
    DEFAULT_SERVERS, ServerConfig, ParallelRunner, Job,
    raise_stack, link_graphs, load_done, parse_phase_timings, build_main,
)

raise_stack()

BIN = "./build/bin/degeneracy_cliques"
TIMEOUT = 1800
OUTCSV = Path("paper_data/bench_local_iterative.csv")
LOGDIR = Path("bench_local_iterative_logs")
FIELDNAMES = ["phase", "graph", "s", "algo", "threads", "status",
              "wall_ms", "total_ms", "peel_ms", "build_ms", "mem_kB", "iter_count"]

# ---- per-server graphs (same split as bench_v3_all.py) ----
SERVER_GRAPHS = {
    "tods1": [
        ("com-dblp",     20),
        ("ca-CondMat",   20),
        ("ca-GrQc",      20),
        ("ca-HepPh",     20),
        ("ca-AstroPh",   20),
        ("com-youtube",  17),
        ("soc-Epinions1",10),
    ],
    "tods2": [
        ("web-Stanford", 40),
        ("dblp-core30",  30),
        ("web-it-2004",  60),
        ("twitter_combined", 15),
        ("web-Google",   25),
    ],
}

# Phase B: parallel-scaling only on the same 3 graphs as paper §7.8
SCALE_GRAPHS = ["com-dblp", "com-youtube", "web-Stanford"]
SCALE_S = {"com-dblp": [3, 5, 10, 15], "com-youtube": [3, 5, 10],
           "web-Stanford": [3, 5, 10, 20, 40]}
SCALE_THREADS = [1, 2, 4, 8, 16, 24]

# Phase A algos
PHASE_A_ALGOS = [
    ("ST_V3",    {"PIVOTER_RUN_ST_V3": "1"}),
    ("LOCAL_V4", {"PIVOTER_RUN_LOCAL_V4": "1"}),
]

# Phase B uses LOCAL_V4 + ST_V3 baseline
PHASE_B_ALGOS = PHASE_A_ALGOS


def gen_phase_a_jobs(graphs: list[tuple[str, int]], done: set):
    """Phase A: T(s) — single thread."""
    for graph, s_max in graphs:
        gpath = f"graphs/{graph}.edges"
        if not os.path.exists(gpath): continue
        for s in range(2, s_max + 1):
            for algo, env in PHASE_A_ALGOS:
                key = ("A", graph, str(s), algo, "1")
                if key in done: continue
                job_env = dict(env)
                job_env["OMP_NUM_THREADS"] = "1"
                yield Job(
                    key=key,
                    cmd=[BIN, gpath, "1", str(s)],
                    env=job_env,
                    log_path=LOGDIR / f"A_{graph}_s{s}_{algo}.log",
                    timeout=TIMEOUT,
                    extra={"phase": "A", "graph": graph, "s": s,
                           "algo": algo, "threads": 1},
                )


def gen_phase_b_jobs(done: set):
    """Phase B: parallel scaling — sequential outer."""
    for graph in SCALE_GRAPHS:
        gpath = f"graphs/{graph}.edges"
        if not os.path.exists(gpath): continue
        for s in SCALE_S[graph]:
            for algo, env in PHASE_B_ALGOS:
                for T in SCALE_THREADS:
                    key = ("B", graph, str(s), algo, str(T))
                    if key in done: continue
                    job_env = dict(env)
                    job_env["OMP_NUM_THREADS"] = str(T)
                    yield Job(
                        key=key,
                        cmd=[BIN, gpath, "1", str(s)],
                        env=job_env,
                        log_path=LOGDIR / f"B_{graph}_s{s}_{algo}_T{T}.log",
                        timeout=TIMEOUT,
                        extra={"phase": "B", "graph": graph, "s": s,
                               "algo": algo, "threads": T},
                    )


def main():
    print("=" * 60)
    print(f"  bench_local_iterative  {time.strftime('%F %T')}")
    print("=" * 60)
    cfg = ServerConfig.detect(DEFAULT_SERVERS)
    print(f"server: {cfg.name}  workers={cfg.max_workers}", flush=True)

    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    # build (only main binary; test_sdct_speed not needed here)
    build_main(["degeneracy_cliques"])

    graphs = SERVER_GRAPHS.get(cfg.name, SERVER_GRAPHS["tods2"])
    avail = link_graphs([g for g, _ in graphs], cfg)
    graphs = [(g, s) for g, s in graphs if g in avail]
    print(f"phase A graphs: {[g for g, _ in graphs]}", flush=True)

    done = load_done(OUTCSV, ("phase", "graph", "s", "algo", "threads"))
    print(f"already done: {len(done)} rows", flush=True)

    runner = ParallelRunner(cfg, OUTCSV, FIELDNAMES)

    def on_finish(job, status, log_text, parsed):
        ex = job.extra
        runner.append_row({
            "phase":   ex["phase"],
            "graph":   ex["graph"],
            "s":       ex["s"],
            "algo":    ex["algo"],
            "threads": ex["threads"],
            "status":  status,
            "wall_ms": f"{parsed['wall_ms']:.1f}",
            "total_ms": f"{parsed['total_ms']:.1f}" if parsed["total_ms"] >= 0 else "",
            "peel_ms":  f"{parsed['peel_ms']:.1f}"  if parsed["peel_ms"]  >= 0 else "",
            "build_ms": f"{parsed['build_ms']:.1f}" if parsed["build_ms"] >= 0 else "",
            "mem_kB":   f"{parsed['mem_kB']:.0f}"   if parsed["mem_kB"]   >= 0 else "",
            "iter_count": parsed["iter_count"] if parsed["iter_count"] >= 0 else "",
        })
        tag = f"[{ex['phase']}] {ex['graph']} s={ex['s']} {ex['algo']} T={ex['threads']}"
        iter_s = f" iter={parsed['iter_count']}" if parsed["iter_count"] > 0 else ""
        peel_s = f" peel={parsed['peel_ms']:.0f}ms" if parsed["peel_ms"] >= 0 else ""
        print(f"  {tag} {status} wall={parsed['wall_ms']:.0f}ms{peel_s}{iter_s}", flush=True)

    # Phase A (parallel outer)
    print("\n=== Phase A: T(s) single-thread ===", flush=True)
    runner.run(gen_phase_a_jobs(graphs, done), on_finish=on_finish)

    # Phase B (sequential outer — one job at a time)
    print("\n=== Phase B: parallel-scaling LOCAL_V4 vs ST_V3 ===", flush=True)
    cfg_seq = ServerConfig(cfg.name, max_workers=1, cpu_load_target=None,
                           mem_limit_gb=cfg.mem_limit_gb,
                           mem_kill_gb=cfg.mem_kill_gb,
                           per_proc_mem_gb=cfg.per_proc_mem_gb,
                           datadir=cfg.datadir)
    runner.cfg = cfg_seq
    done = load_done(OUTCSV, ("phase", "graph", "s", "algo", "threads"))
    runner.run(gen_phase_b_jobs(done), on_finish=on_finish)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
