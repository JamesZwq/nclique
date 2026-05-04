#!/usr/bin/env python3
"""
LocalH (LOCAL_V4) single-thread experiment for VLDB R=1 paper §7 RQ5.

LocalH is the synchronous fixed-point solver — direct h-index iteration over
the implicit hypergraph H_s(G).  It is the semantic baseline for PeelH (V3),
not a performance contender: paper claims LOCAL_V4 takes 396 s on com-dblp
at s=10 vs 7 ms for V3 peel.  This script collects the raw numbers behind
that claim across a representative subset of graphs.

Each invocation runs `degeneracy_cliques` (the standard binary) with
PIVOTER_RUN_LOCAL_V4=1, single-threaded, wrapped in /usr/bin/time -v so
peak RSS, page faults, and exit status land in the log.

Usage:
    nohup python3 bench_local_v4.py tods2 > /tmp/bench_local_v4.log 2>&1 &
"""
from __future__ import annotations
import csv as _csv
import os, re, subprocess, sys, time
from pathlib import Path

from bench_lib import (
    DEFAULT_SERVERS, ServerConfig, raise_stack, link_graphs, build_main,
)

raise_stack()

BIN        = "./build/bin/degeneracy_cliques"
TIME_BIN   = "/usr/bin/time"
TIMEOUT    = 1800            # per-cell upper bound (30 min)
OUTCSV     = Path("paper_data/bench_local_v4.csv")
LOGDIR     = Path("bench_local_v4_logs")
RUNS_PER_CFG = 1             # LocalH is deterministic at T=1; no noise floor

FIELDNAMES = [
    "graph", "s", "algo", "status", "wall_ms", "took_ms",
    "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
    "time_pagefaults_major", "time_pagefaults_minor",
    "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status",
]

# Representative subset: paper RQ5 wants the LocalH/PeelH gap *across the
# regime*, not exhaustive coverage. Skip web-it-2004 (T=1 sequential SDCT
# already takes 90 s; LocalH on top would push past TIMEOUT for every s).
SERVER_GRAPHS = {
    "tods2": [
        ("com-amazon.ungraph",       [3, 5, 8]),
        ("com-dblp",                 [3, 5, 8, 10]),     # s=10 = canonical paper number
        ("com-youtube",              [3, 5, 8]),
        ("soc-pokec-relationships",  [3, 5]),
        ("twitter_combined",         [3]),
        ("wiki-Talk",                [3, 5]),
        ("web-Google",               [3, 5, 8]),
        ("web-Stanford",             [3, 5, 10]),
        ("dblp-core30",              [3, 5, 10]),
    ],
}

_NUM     = r"[\d.]+(?:e[+-]?\d+)?"
_TOOK_RE = re.compile(rf"Local H-index V4 r=1 took:\s+({_NUM})\s+ms")


def parse_took(stdout: str) -> float | None:
    for line in stdout.splitlines():
        m = _TOOK_RE.search(line)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                return None
    return None


def parse_time_v(stderr: str) -> dict:
    out = {k: "" for k in (
        "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
        "time_pagefaults_major", "time_pagefaults_minor",
        "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status")}
    for raw in stderr.splitlines():
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


def load_existing_keys(csv_path: Path) -> set:
    """Skip rows already at OK status — LocalH at T=1 is deterministic."""
    if not csv_path.exists(): return set()
    done = set()
    with csv_path.open() as f:
        for row in _csv.DictReader(f):
            if row.get("status") == "OK":
                try:
                    done.add((row["graph"], int(row["s"])))
                except (ValueError, KeyError):
                    pass
    return done


def append_row(csv_path: Path, row: dict) -> None:
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("a", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=FIELDNAMES)
        if write_header: w.writeheader()
        w.writerow({k: row.get(k, "") for k in FIELDNAMES})


def run_one(graph: str, s: int, run_idx: int) -> dict:
    gpath = f"graphs/{graph}.edges"
    log_path = LOGDIR / f"{graph}_s{s}_r{run_idx}.log"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["PIVOTER_RUN_LOCAL_V4"] = "1"
    args = [TIME_BIN, "-v", BIN, gpath, "1", str(s)]
    cmd_line = " ".join(args)
    t0 = time.time()
    try:
        proc = subprocess.run(args, env=env, capture_output=True, text=True,
                              timeout=TIMEOUT)
        wall_ms = (time.time() - t0) * 1000.0
        log_path.write_text(
            f"[CMD] {cmd_line}\n"
            f"[ENV] OMP_NUM_THREADS=1 PIVOTER_RUN_LOCAL_V4=1\n"
            f"[RC]  {proc.returncode}\n\n"
            f"--STDOUT--\n{proc.stdout}\n"
            f"--STDERR--\n{proc.stderr}\n"
        )
        time_v = parse_time_v(proc.stderr)
        if proc.returncode == 0:
            took = parse_took(proc.stdout)
            return {"status": "OK", "wall_ms": wall_ms,
                    "took_ms": (took if took is not None else ""),
                    **time_v}
        else:
            return {"status": f"ERROR({proc.returncode})",
                    "wall_ms": wall_ms, **time_v}
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode("utf-8", "replace") if e.stdout else ""
        err = e.stderr.decode("utf-8", "replace") if e.stderr else ""
        log_path.write_text(
            f"[CMD] {cmd_line}\n"
            f"[ENV] OMP_NUM_THREADS=1 PIVOTER_RUN_LOCAL_V4=1\n"
            f"[TIMEOUT after {TIMEOUT}s]\n\n"
            f"--STDOUT--\n{out}\n--STDERR--\n{err}\n"
        )
        return {"status": "TIMEOUT", "wall_ms": TIMEOUT * 1000.0,
                **parse_time_v(err)}


def main():
    print("=" * 60)
    print(f"  bench_local_v4  {time.strftime('%F %T')}")
    print("=" * 60)
    cfg = ServerConfig.detect(DEFAULT_SERVERS)
    print(f"server: {cfg.name}", flush=True)

    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    build_main(["degeneracy_cliques"])

    graphs = SERVER_GRAPHS.get(cfg.name, SERVER_GRAPHS["tods2"])
    avail = link_graphs([g for g, _ in graphs], cfg)
    graphs = [(g, ss) for g, ss in graphs if g in avail]
    print(f"graphs: {[g for g, _ in graphs]}", flush=True)

    skip = load_existing_keys(OUTCSV)
    print(f"skipping {len(skip)} fully-done (graph,s) cells", flush=True)

    total = sum(len(ss) for _, ss in graphs) * RUNS_PER_CFG
    done_cnt = 0
    for graph, s_list in graphs:
        for s in s_list:
            if (graph, s) in skip:
                done_cnt += RUNS_PER_CFG
                continue
            for run_idx in range(RUNS_PER_CFG):
                r = run_one(graph, s, run_idx)
                row = {
                    "graph": graph, "s": s, "algo": "LOCAL_V4",
                    "status":  r["status"],
                    "wall_ms": f"{r['wall_ms']:.1f}",
                    "took_ms": r.get("took_ms", ""),
                    **{k: r.get(k, "") for k in (
                        "time_max_rss_kB","time_user_sec","time_sys_sec",
                        "time_elapsed","time_pagefaults_major",
                        "time_pagefaults_minor","time_voluntary_ctxt",
                        "time_involuntary_ctxt","time_exit_status")},
                }
                append_row(OUTCSV, row)
                done_cnt += 1
                print(f"  [{done_cnt}/{total}] {graph} s={s} run={run_idx} "
                      f"{r['status']} took={r.get('took_ms', 'N/A')}ms "
                      f"wall={r['wall_ms']/1000:.1f}s", flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
