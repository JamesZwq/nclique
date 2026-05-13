#!/usr/bin/env python3
"""
(1,s)-nucleus hierarchy bench: time + memory across all 10 paper graphs.

Compares two algos using PIVOTER_RUN_ST_V3=1 (the SPIN-star pipeline):
  - ST_V3      : no hierarchy build (baseline)
  - ST_V3_HIER : PIVOTER_DUMP_HIER=<sink> triggers post-peel BuildHier

For every (graph, s, algo) we record
  graph, s, algo, status, wall_ms, total_ms, build_ms, peel_ms, hier_ms, mem_kB

Graphs are run small-to-large by edge count.  Workers cap at 32 per the
user instruction; the ParallelRunner mem-eviction kills the newest job
when system memory exceeds the configured budget.

Usage:
    nohup python3 bench_r1_hierarchy.py tods1 > bench_r1_hierarchy.log 2>&1 &
"""
from __future__ import annotations

import csv
import os
import re
import sys
import tempfile
from pathlib import Path

from bench_lib import (
    DEFAULT_SERVERS, ServerConfig, ParallelRunner, Job,
    raise_stack, link_graphs, load_done, build_main, parse_phase_timings,
)

raise_stack()

BIN     = "./build/bin/degeneracy_cliques"
TIMEOUT = 1800
OUTCSV  = Path("paper_data/bench_r1_hierarchy.csv")
LOGDIR  = Path("bench_r1_hierarchy_logs")
FIELDNAMES = ["graph", "s", "algo", "status",
              "wall_ms", "total_ms", "build_ms", "peel_ms",
              "hier_ms", "mem_kB"]

# Smallest -> largest by edge count (matches paper Section 8 dataset table).
GRAPHS = [
    "com-amazon.ungraph",       # 0.93M edges
    "com-dblp",                 # 1.05M edges
    "twitter_combined",         # 1.34M edges
    "web-Stanford",             # 1.99M edges
    "com-youtube",              # 2.99M edges
    "web-Google",               # 4.32M edges
    "wiki-Talk",                # 4.66M edges
    "web-it-2004",              # 7.18M edges
    "soc-pokec-relationships",  # 22.3M edges
    "com-orkut",                # 117M edges
]

# Per-graph s cap (matches bench_r1_main.py).  We sweep the standard set
# {3,5,7,9,11,13,15} and skip any value above the per-graph cap.
S_VALUES = [3, 5, 7, 9, 11, 13, 15]
S_MAX = {
    "com-dblp": 30, "com-amazon.ungraph": 30, "twitter_combined": 30,
    "web-Stanford": 61, "web-Google": 40, "com-youtube": 17,
    "web-it-2004": 430, "wiki-Talk": 25, "soc-pokec-relationships": 25,
    "com-orkut": 15,
}

# Four algos.  HIER_SINK is a process-local sink file (per-job, not read
# back); we only measure the emit time + memory.
#
#   ST_V3        : SPIN-star baseline, no hierarchy.
#   ST_V3_HIER   : SPIN-star + post-peel BuildHier (our O(Sigma log Sigma)
#                  one-pass routine).
#   REF_R1       : NuclearCD baseline (mutable-CPI), no hierarchy.
#   REF_R1_HIER  : NuclearCD baseline + CND-style level-DSU hierarchy
#                  (paper baseline for Section 7).
ALGOS = [
    ("ST_V3",       {"PIVOTER_RUN_ST_V3": "1"}),
    ("ST_V3_HIER",  {"PIVOTER_RUN_ST_V3": "1",
                     "PIVOTER_DUMP_HIER":     "PER_JOB_SINK"}),
    ("REF_R1",      {}),
    ("REF_R1_HIER", {"PIVOTER_DUMP_HIER_REF": "PER_JOB_SINK"}),
]


_RE_HIER     = re.compile(r"hier r=1[^:]*took:\s*([\d.]+)\s*ms")
_RE_HIER_REF = re.compile(r"hier r=1 \(CND level-DSU\) took:\s*([\d.]+)\s*ms")


def parse_extra(txt: str) -> dict:
    out = parse_phase_timings(txt)
    # Prefer the CND-style match if present (its regex is more specific);
    # otherwise fall back to the generic ST_V3 hier line.
    m = _RE_HIER_REF.search(txt) or _RE_HIER.search(txt)
    out["hier_ms"] = float(m.group(1)) if m else -1.0
    return out


def gen_jobs(graphs: list[str], done: set):
    for graph in graphs:
        gpath = f"graphs/{graph}.edges"
        if not os.path.exists(gpath):
            print(f"[skip] missing {gpath}", flush=True)
            continue
        smax = S_MAX.get(graph, 15)
        for s in S_VALUES:
            if s > smax:
                continue
            for algo, env in ALGOS:
                key = (graph, str(s), algo)
                if key in done:
                    continue
                # Per-job sink so concurrent jobs don't trample one file.
                env_copy = dict(env)
                for ek in ("PIVOTER_DUMP_HIER", "PIVOTER_DUMP_HIER_REF"):
                    if env_copy.get(ek) == "PER_JOB_SINK":
                        env_copy[ek] = tempfile.mktemp(
                            prefix=f"hier_{graph}_s{s}_{algo}_", suffix=".tsv")
                env_copy["OMP_NUM_THREADS"] = "1"
                yield Job(
                    key=key,
                    cmd=[BIN, gpath, "1", str(s), "degen"],
                    env=env_copy,
                    log_path=LOGDIR / f"{graph}_s{s}_{algo}.log",
                    timeout=TIMEOUT,
                )


def on_finish(job: Job, status: str, log_text: str, parsed: dict):
    graph, s, algo = job.key
    parsed.update(parse_extra(log_text))
    row = {
        "graph": graph, "s": s, "algo": algo, "status": status,
        "wall_ms":  f"{parsed.get('wall_ms', -1.0):.1f}",
        "total_ms": f"{parsed.get('total_ms', -1.0):.1f}",
        "build_ms": f"{parsed.get('build_ms', -1.0):.1f}",
        "peel_ms":  f"{parsed.get('peel_ms', -1.0):.1f}",
        "hier_ms":  f"{parsed.get('hier_ms', -1.0):.1f}",
        "mem_kB":   f"{parsed.get('mem_kB', -1.0):.0f}",
    }
    runner.append_row(row)
    # Clean up per-job sink file if any (we never read it back).
    for ek in ("PIVOTER_DUMP_HIER", "PIVOTER_DUMP_HIER_REF"):
        sink = job.env.get(ek)
        if sink and sink not in ("/dev/null", "PER_JOB_SINK"):
            try: os.unlink(sink)
            except OSError: pass


def announce(job: Job, status: str, parsed: dict):
    g, s, algo = job.key
    print(f"  [done] {g:24s} s={s:>2}  {algo:11s}  "
          f"status={status}  total={parsed.get('total_ms', -1):.0f}ms  "
          f"peel={parsed.get('peel_ms', -1):.0f}ms  "
          f"hier={parsed.get('hier_ms', -1):.1f}ms  "
          f"mem={parsed.get('mem_kB', -1)/1024:.0f}MB",
          flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in DEFAULT_SERVERS:
        print(f"Usage: {sys.argv[0]} <{ '|'.join(DEFAULT_SERVERS) }>", file=sys.stderr)
        sys.exit(1)
    cfg = DEFAULT_SERVERS[sys.argv[1]]
    # User-requested 32-worker cap.
    cfg = ServerConfig(cfg.name, max_workers=min(cfg.max_workers, 32),
                       cpu_load_target=cfg.cpu_load_target,
                       mem_limit_gb=cfg.mem_limit_gb,
                       mem_kill_gb=cfg.mem_kill_gb,
                       per_proc_mem_gb=cfg.per_proc_mem_gb)
    print(f"[cfg] {cfg.name} workers={cfg.max_workers} mem_limit={cfg.mem_limit_gb}GB "
          f"mem_kill={cfg.mem_kill_gb}GB", flush=True)

    LOGDIR.mkdir(parents=True, exist_ok=True)
    available = link_graphs(GRAPHS, cfg)
    print(f"[graphs] {len(available)}/{len(GRAPHS)} available", flush=True)

    build_main(["degeneracy_cliques"])

    done = load_done(OUTCSV, key_fields=("graph", "s", "algo"))
    runner = ParallelRunner(cfg, OUTCSV, FIELDNAMES)
    jobs = gen_jobs(available, done)

    print(f"[start] writing {OUTCSV}", flush=True)
    runner.run(jobs, on_finish=on_finish, announce_done=announce)
    print("[done]", flush=True)
