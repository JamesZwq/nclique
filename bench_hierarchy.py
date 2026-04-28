#!/usr/bin/env python3
"""
Experiment 5 — Hierarchy emission timing (Theorem 4: hierarchy as peel
by-product).

For r>=3 the hierarchy is already exposed via PIVOTER_RUN_REGION_V3LM_HIER,
which extends V3LM with class-based DSU post-processing. We measure the
DSU overhead by comparing V3LM (no hierarchy) vs V3LM_HIER (with) — the
paper's claim is that the post-processing is ~0-1% of peel.

For r=1 the equivalent flag is not yet wired into the centralised peel.
Section §6.1 promises Theorem 4 holds at r=1 by the same argument as r>=3
but the overhead has not been measured. To collect r=1 data we need a
small change in NCliqueVertexCoreDecompositionST_V3.cpp — emit (parent,
child) edges of the hierarchy DAG when PIVOTER_DUMP_HIERARCHY=1 is set.
That extension is tracked separately; this script measures the existing
r>=3 hierarchy overhead.

Usage:
    nohup python3 bench_hierarchy.py tods2 > bench_hierarchy.log 2>&1 &
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
OUTCSV = Path("paper_data/bench_hierarchy.csv")
LOGDIR = Path("bench_hierarchy_logs")
FIELDNAMES = ["graph", "r", "s", "algo", "status",
              "wall_ms", "total_ms", "peel_ms", "hier_ms", "mem_kB"]

# r=3 only here; if you have probed max-clique sizes per graph, you can
# extend the (r,s) list per graph below to cover r=4, r=5, ...
SERVER_GRAPHS = {
    "tods1": [
        ("com-dblp",     [(3, 5), (3, 8), (4, 5), (4, 8)]),
        ("com-youtube",  [(3, 5), (4, 5), (4, 8)]),
        ("ca-CondMat",   [(3, 5), (4, 5), (5, 6)]),
    ],
    "tods2": [
        ("web-Stanford", [(3, 5), (3, 8), (4, 5)]),
        ("dblp-core30",  [(3, 5), (4, 5), (5, 6)]),
        ("web-it-2004",  [(3, 5), (3, 8)]),
    ],
}

ALGOS = [
    ("V3LM",      {"PIVOTER_RUN_REGION_V3LM":      "1"}),
    ("V3LM_HIER", {"PIVOTER_RUN_REGION_V3LM_HIER": "1"}),
]


def parse_hier(txt: str) -> float:
    import re
    m = re.search(r'Hierarchy post-processing(?: \(class-based\))?:\s*([\d.]+)', txt)
    return float(m.group(1)) if m else -1.0


def gen_jobs(graphs, done):
    for graph, rs_list in graphs:
        gpath = f"graphs/{graph}.edges"
        if not os.path.exists(gpath): continue
        for r, s in rs_list:
            for algo, env in ALGOS:
                key = (graph, str(r), str(s), algo)
                if key in done: continue
                yield Job(
                    key=key,
                    cmd=[BIN, gpath, str(r), str(s), "degen"],
                    env=dict(env, OMP_NUM_THREADS="1"),
                    log_path=LOGDIR / f"{graph}_r{r}_s{s}_{algo}.log",
                    timeout=TIMEOUT,
                    extra={"graph": graph, "r": r, "s": s, "algo": algo},
                )


def main():
    print("=" * 60)
    print(f"  bench_hierarchy  {time.strftime('%F %T')}")
    print("=" * 60)
    cfg = ServerConfig.detect(DEFAULT_SERVERS)
    print(f"server: {cfg.name}  workers={cfg.max_workers}", flush=True)

    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    build_main(["degeneracy_cliques"])

    graphs = SERVER_GRAPHS.get(cfg.name, SERVER_GRAPHS["tods2"])
    avail = link_graphs([g for g, _ in graphs], cfg)
    graphs = [(g, rs) for g, rs in graphs if g in avail]
    print(f"graphs: {[g for g, _ in graphs]}", flush=True)

    done = load_done(OUTCSV, ("graph", "r", "s", "algo"))
    print(f"already done: {len(done)} rows", flush=True)

    runner = ParallelRunner(cfg, OUTCSV, FIELDNAMES)

    def on_finish(job, status, log_text, parsed):
        ex = job.extra
        hier_ms = parse_hier(log_text)
        runner.append_row({
            "graph":   ex["graph"], "r": ex["r"], "s": ex["s"], "algo": ex["algo"],
            "status":  status,
            "wall_ms": f"{parsed['wall_ms']:.1f}",
            "total_ms": f"{parsed['total_ms']:.1f}" if parsed["total_ms"] >= 0 else "",
            "peel_ms":  f"{parsed['peel_ms']:.1f}"  if parsed["peel_ms"]  >= 0 else "",
            "hier_ms":  f"{hier_ms:.1f}"            if hier_ms             >= 0 else "",
            "mem_kB":   f"{parsed['mem_kB']:.0f}"   if parsed["mem_kB"]   >= 0 else "",
        })
        tag = f"{ex['graph']:>14} r={ex['r']} s={ex['s']:>3} {ex['algo']:>10}"
        h_str = f" hier={hier_ms:.1f}ms" if hier_ms >= 0 else ""
        print(f"  {tag} {status} wall={parsed['wall_ms']:.0f}ms{h_str}", flush=True)

    runner.run(gen_jobs(graphs, done), on_finish=on_finish)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
