#!/usr/bin/env python3
"""
SPIN (LocalH / LOCAL_V4) full sweep — matches the v3 main bench coverage.

Differences from bench_local_v4.py:
  * All 10 graphs at full s ranges (not the curated subset).
  * 1h timeout (was 30 min).
  * Skip-floor: once a graph times out at some s, skip higher s for that graph.
  * Multi-threaded via ParallelRunner with memory gate (was serial).

Resume-friendly: skips (graph, s) cells already at OK status.

Usage:
    nohup python3 bench_local_v4_full.py tods2 > /tmp/bench_local_v4_full.log 2>&1 &
"""
from __future__ import annotations
import csv as _csv, os, time
from pathlib import Path
from bench_lib import (
    DEFAULT_SERVERS, ServerConfig, raise_stack, link_graphs, build_main,
    ParallelRunner, Job, load_done,
)

raise_stack()

BIN      = "./build/bin/degeneracy_cliques"
TIME_BIN = "/usr/bin/time"
TIMEOUT  = 3600
OUTCSV   = Path("paper_data/bench_local_v4.csv")
LOGDIR   = Path("bench_local_v4_logs")

FIELDNAMES = [
    "graph", "s", "algo", "status", "wall_ms", "took_ms",
    "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
    "time_pagefaults_major", "time_pagefaults_minor",
    "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status",
]

# Match v3 main bench coverage so SPIN★ vs SPIN comparison is comparable.
GRAPHS: list[tuple[str, int]] = [
    ("com-dblp",                 30),
    ("com-amazon.ungraph",       30),
    ("twitter_combined",         30),
    ("web-Stanford",             61),
    ("web-Google",               40),
    ("com-youtube",              17),
    ("web-it-2004",             430),
    ("wiki-Talk",                25),
    ("soc-pokec-relationships",  25),
    ("com-orkut",                15),
    ("dblp-core30",              30),
]


def main():
    print("=" * 60)
    print(f"  bench_local_v4_full  {time.strftime('%F %T')}")
    print("=" * 60)
    cfg = ServerConfig.detect(DEFAULT_SERVERS)
    print(f"server: {cfg.name}  max_workers={cfg.max_workers}  "
          f"mem_limit={cfg.mem_limit_gb}GB", flush=True)

    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)
    build_main(["degeneracy_cliques"])

    avail = link_graphs([g for g, _ in GRAPHS], cfg)
    graphs = [(g, smax) for g, smax in GRAPHS if g in avail]
    print(f"graphs: {[g for g, _ in graphs]}", flush=True)

    # Existing OK cells — skip
    done_ok = load_done(OUTCSV, ("graph", "s"))
    # Existing TIMEOUT/OOM cells — also skip and trigger skip-floor for that graph
    done_skip = set()
    skip_floor: dict[str, int] = {}
    if OUTCSV.exists():
        with OUTCSV.open() as f:
            for row in _csv.DictReader(f):
                try:
                    g, s = row["graph"], int(row["s"])
                except (KeyError, ValueError):
                    continue
                st = row.get("status", "")
                if st in ("TIMEOUT", "OOM") or st.startswith("ERROR"):
                    done_skip.add((g, str(s)))
                    cur = skip_floor.get(g)
                    if cur is None or s < cur:
                        skip_floor[g] = s
    print(f"already OK: {len(done_ok)}; prior TIMEOUT/OOM: {len(done_skip)}", flush=True)
    print(f"skip-floor (graph -> first failed s): {skip_floor}", flush=True)

    # Build job list
    jobs: list[Job] = []
    for g, smax in graphs:
        for s in range(2, smax + 1):
            if (g, str(s)) in done_ok or (g, str(s)) in done_skip:
                continue
            sf = skip_floor.get(g)
            if sf is not None and s >= sf:
                continue
            log_path = LOGDIR / f"{g}_s{s}.log"
            jobs.append(Job(
                key=(g, s),
                cmd=[TIME_BIN, "-v", BIN, f"graphs/{g}.edges", "1", str(s)],
                env={"OMP_NUM_THREADS": "1", "PIVOTER_RUN_LOCAL_V4": "1"},
                log_path=log_path,
                timeout=TIMEOUT,
            ))
    print(f"[plan] {len(jobs)} jobs to schedule", flush=True)

    # Skip-floor enforcement during run: if a job for graph G times out at s,
    # all higher-s jobs for G already in the queue should be skipped. Easiest
    # is to prune the iterable just-in-time.
    def gen_jobs():
        for job in jobs:
            g, s = job.key
            sf = skip_floor.get(g)
            if sf is not None and s >= sf:
                continue
            yield job

    runner = ParallelRunner(cfg, OUTCSV, FIELDNAMES)
    done_cnt = [0]
    total = len(jobs)

    def parse_took(log_text: str) -> str:
        import re
        m = re.search(r"Local H-index V4 r=1 took:\s+([\d.]+)\s+ms", log_text)
        return m.group(1) if m else ""

    def parse_time_v(log_text: str) -> dict:
        out = {k: "" for k in (
            "time_max_rss_kB","time_user_sec","time_sys_sec","time_elapsed",
            "time_pagefaults_major","time_pagefaults_minor",
            "time_voluntary_ctxt","time_involuntary_ctxt","time_exit_status")}
        for raw in log_text.splitlines():
            line = raw.strip()
            if   line.startswith("Maximum resident set size"):
                out["time_max_rss_kB"]       = line.split(":")[-1].strip()
            elif line.startswith("User time (seconds)"):
                out["time_user_sec"]         = line.split(":")[-1].strip()
            elif line.startswith("System time (seconds)"):
                out["time_sys_sec"]          = line.split(":")[-1].strip()
            elif line.startswith("Elapsed (wall clock) time"):
                out["time_elapsed"]          = line.split(": ", 1)[-1].strip()
            elif line.startswith("Major (requiring I/O) page faults"):
                out["time_pagefaults_major"] = line.split(":")[-1].strip()
            elif line.startswith("Minor (reclaiming a frame) page faults"):
                out["time_pagefaults_minor"] = line.split(":")[-1].strip()
            elif line.startswith("Voluntary context switches"):
                out["time_voluntary_ctxt"]   = line.split(":")[-1].strip()
            elif line.startswith("Involuntary context switches"):
                out["time_involuntary_ctxt"] = line.split(":")[-1].strip()
            elif line.startswith("Exit status"):
                out["time_exit_status"]      = line.split(":")[-1].strip()
        return out

    def on_finish(job: Job, status: str, log_text: str, parsed: dict):
        g, s = job.key
        row = {
            "graph": g, "s": s, "algo": "LOCAL_V4", "status": status,
            "wall_ms": f"{parsed.get('wall_ms', 0):.1f}",
            "took_ms": parse_took(log_text),
            **parse_time_v(log_text),
        }
        runner.append_row(row)
        done_cnt[0] += 1
        print(f"  [{done_cnt[0]}/{total}] {g} s={s} {status} "
              f"took={row['took_ms']}ms wall={parsed.get('wall_ms',0)/1000:.1f}s "
              f"rss={row['time_max_rss_kB']}kB", flush=True)
        # skip-floor on this graph
        if status in ("TIMEOUT", "OOM") or status.startswith("ERROR"):
            cur = skip_floor.get(g)
            if cur is None or s < cur:
                skip_floor[g] = s
                print(f"    -> skip-floor: {g} s>={s}", flush=True)

    runner.run(gen_jobs(), on_finish)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
