#!/usr/bin/env python3
"""
Experiment 4 — Vertex-ordering ablation (Experimental.tex auxiliary).

Sweeps the binary's <sort_option> argument across degen / degenR /
degree / degreeR / default for ST_V3 (the optimised peel) on the same
graphs as the ablation experiment. Records build / peel / mem so the
paper can quantify how much the degeneracy ordering buys vs simpler
alternatives.

Usage:
    nohup python3 bench_sort_options.py tods2 > bench_sort_options.log 2>&1 &
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
OUTCSV = Path("paper_data/bench_sort_options.csv")
LOGDIR = Path("bench_sort_options_logs")
FIELDNAMES = ["graph", "s", "sort", "status",
              "wall_ms", "total_ms", "peel_ms", "build_ms", "mem_kB"]

SERVER_GRAPHS = {
    "tods1": [
        ("com-dblp",     [3, 5, 8, 10, 15]),
        ("com-youtube",  [3, 5, 8, 12]),
        ("ca-CondMat",   [3, 5, 10, 15]),
    ],
    "tods2": [
        ("web-Stanford", [3, 5, 10, 20, 40]),
        ("dblp-core30",  [3, 5, 10, 20]),
        ("web-it-2004",  [3, 10, 30, 100]),
    ],
}

SORTS = ["degen", "degenR", "degree", "degreeR", "default"]


def gen_jobs(graphs, done):
    for graph, s_list in graphs:
        gpath = f"graphs/{graph}.edges"
        if not os.path.exists(gpath): continue
        for s in s_list:
            for sort in SORTS:
                key = (graph, str(s), sort)
                if key in done: continue
                yield Job(
                    key=key,
                    cmd=[BIN, gpath, "1", str(s), sort],
                    env={"PIVOTER_RUN_ST_V3": "1", "OMP_NUM_THREADS": "1"},
                    log_path=LOGDIR / f"{graph}_s{s}_{sort}.log",
                    timeout=TIMEOUT,
                    extra={"graph": graph, "s": s, "sort": sort},
                )


def main():
    print("=" * 60)
    print(f"  bench_sort_options  {time.strftime('%F %T')}")
    print("=" * 60)
    cfg = ServerConfig.detect(DEFAULT_SERVERS)
    print(f"server: {cfg.name}  workers={cfg.max_workers}", flush=True)

    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    build_main(["degeneracy_cliques"])

    graphs = SERVER_GRAPHS.get(cfg.name, SERVER_GRAPHS["tods2"])
    avail = link_graphs([g for g, _ in graphs], cfg)
    graphs = [(g, ss) for g, ss in graphs if g in avail]
    print(f"graphs: {[g for g, _ in graphs]}", flush=True)

    done = load_done(OUTCSV, ("graph", "s", "sort"))
    print(f"already done: {len(done)} rows", flush=True)

    runner = ParallelRunner(cfg, OUTCSV, FIELDNAMES)

    def on_finish(job, status, log_text, parsed):
        ex = job.extra
        runner.append_row({
            "graph":   ex["graph"], "s": ex["s"], "sort": ex["sort"],
            "status":  status,
            "wall_ms": f"{parsed['wall_ms']:.1f}",
            "total_ms": f"{parsed['total_ms']:.1f}" if parsed["total_ms"] >= 0 else "",
            "peel_ms":  f"{parsed['peel_ms']:.1f}"  if parsed["peel_ms"]  >= 0 else "",
            "build_ms": f"{parsed['build_ms']:.1f}" if parsed["build_ms"] >= 0 else "",
            "mem_kB":   f"{parsed['mem_kB']:.0f}"   if parsed["mem_kB"]   >= 0 else "",
        })
        tag = f"{ex['graph']:>14} s={ex['s']:>3} {ex['sort']:>8}"
        peel_s = f" peel={parsed['peel_ms']:.1f}ms" if parsed["peel_ms"] >= 0 else ""
        print(f"  {tag} {status} wall={parsed['wall_ms']:.0f}ms{peel_s}", flush=True)

    runner.run(gen_jobs(graphs, done), on_finish=on_finish)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
