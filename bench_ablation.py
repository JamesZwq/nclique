#!/usr/bin/env python3
"""
Experiment 2 — Controlled ablation ladder (Experimental.tex §6.5).

Four-step incremental ladder; each step adds one optimisation:

    REF   : mutable-tree baseline (no env flags)
    ST    : + immutable tree + nCr-delta + flat CSR
    ST_V2 : + tree-free build (PIVOTER_RUN_ST_V2)
    ST_V3 : + sparse bucket queue (PIVOTER_RUN_ST_V3)

Per-graph speedup of ST/V2/V3 over REF is reported on the same paper-3:
    com-dblp, com-youtube, web-Stanford
plus extras filling out the per-server pool.

Usage:
    nohup python3 bench_ablation.py tods2 > bench_ablation.log 2>&1 &
"""
from __future__ import annotations
import os, time
from pathlib import Path

from bench_lib import (
    DEFAULT_SERVERS, ServerConfig, ParallelRunner, Job,
    raise_stack, link_graphs, load_done, build_main,
)

raise_stack()

BIN = "./build/bin/degeneracy_cliques"
TIMEOUT = 1800
OUTCSV = Path("paper_data/bench_ablation.csv")
LOGDIR = Path("bench_ablation_logs")
FIELDNAMES = ["graph", "s", "algo", "status",
              "wall_ms", "total_ms", "peel_ms", "build_ms", "mem_kB"]

SERVER_GRAPHS = {
    "tods1": [
        ("com-dblp",     30),
        ("com-youtube",  17),
        ("ca-CondMat",   20),
        ("ca-HepPh",     22),
        ("ca-AstroPh",   25),
    ],
    "tods2": [
        ("web-Stanford", 60),
        ("dblp-core30",  30),
        ("web-Google",   25),
        ("web-it-2004",  60),
    ],
}

ALGOS = [
    ("REF",   {}),
    ("ST",    {"PIVOTER_RUN_ST":    "1"}),
    ("ST_V2", {"PIVOTER_RUN_ST_V2": "1"}),
    ("ST_V3", {"PIVOTER_RUN_ST_V3": "1"}),
]


def gen_jobs(graphs: list[tuple[str, int]], done: set):
    for graph, s_max in graphs:
        gpath = f"graphs/{graph}.edges"
        if not os.path.exists(gpath): continue
        for s in range(2, s_max + 1):
            for algo, env in ALGOS:
                key = (graph, str(s), algo)
                if key in done: continue
                yield Job(
                    key=key,
                    cmd=[BIN, gpath, "1", str(s)],
                    env=dict(env, OMP_NUM_THREADS="1"),
                    log_path=LOGDIR / f"{graph}_s{s}_{algo}.log",
                    timeout=TIMEOUT,
                    extra={"graph": graph, "s": s, "algo": algo},
                )


def main():
    print("=" * 60)
    print(f"  bench_ablation  {time.strftime('%F %T')}")
    print("=" * 60)
    cfg = ServerConfig.detect(DEFAULT_SERVERS)
    print(f"server: {cfg.name}  workers={cfg.max_workers}", flush=True)

    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    build_main(["degeneracy_cliques"])

    graphs = SERVER_GRAPHS.get(cfg.name, SERVER_GRAPHS["tods2"])
    avail = link_graphs([g for g, _ in graphs], cfg)
    graphs = [(g, s) for g, s in graphs if g in avail]
    print(f"graphs: {[g for g, _ in graphs]}", flush=True)

    done = load_done(OUTCSV, ("graph", "s", "algo"))
    print(f"already done: {len(done)} rows", flush=True)

    runner = ParallelRunner(cfg, OUTCSV, FIELDNAMES)

    def on_finish(job, status, log_text, parsed):
        ex = job.extra
        runner.append_row({
            "graph":   ex["graph"], "s": ex["s"], "algo": ex["algo"],
            "status":  status,
            "wall_ms": f"{parsed['wall_ms']:.1f}",
            "total_ms": f"{parsed['total_ms']:.1f}" if parsed["total_ms"] >= 0 else "",
            "peel_ms":  f"{parsed['peel_ms']:.1f}"  if parsed["peel_ms"]  >= 0 else "",
            "build_ms": f"{parsed['build_ms']:.1f}" if parsed["build_ms"] >= 0 else "",
            "mem_kB":   f"{parsed['mem_kB']:.0f}"   if parsed["mem_kB"]   >= 0 else "",
        })
        tag = f"{ex['graph']:>14} s={ex['s']:>3} {ex['algo']:>5}"
        peel_s = f" peel={parsed['peel_ms']:.1f}ms" if parsed["peel_ms"] >= 0 else ""
        print(f"  {tag} {status} wall={parsed['wall_ms']:.0f}ms{peel_s}", flush=True)

    runner.run(gen_jobs(graphs, done), on_finish=on_finish)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
