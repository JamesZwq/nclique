#!/usr/bin/env python3
"""
Experiment 3 — Parallel SDCT walk thread-scaling (Experimental.tex §6.6,
Cor 1).

Times SDCT_Par4 (winner from §6.6 ablation) at OMP_NUM_THREADS in
{1,2,4,8,16,24} on com-dblp, com-youtube, web-Stanford for the same s
range used in the breakdown experiment.

Outer loop is sequential (one job at a time) so each binary call binds all
cores. Build of the auxiliary binary `bench_sdct_par_scale` is checked
first; that binary parses (graph, s, threads) and emits one PAR4 line.

Usage:
    nohup python3 bench_par_sdct.py tods2 > bench_par_sdct.log 2>&1 &
"""
from __future__ import annotations
import os, re, subprocess, sys, time
from pathlib import Path

from bench_lib import (
    DEFAULT_SERVERS, ServerConfig, raise_stack, link_graphs, build_main,
)

raise_stack()

BIN = "./build/bin/bench_sdct_par_scale"
TIMEOUT = 1800
OUTCSV = Path("paper_data/bench_par_sdct.csv")
LOGDIR = Path("bench_par_sdct_logs")
FIELDNAMES = ["graph", "s", "threads", "status", "wall_ms", "build_ms",
              "leaves", "ref_leaves", "verify",
              # /usr/bin/time -v fields, parsed from stderr
              "time_max_rss_kB", "time_user_sec", "time_sys_sec",
              "time_elapsed", "time_pagefaults_major", "time_pagefaults_minor",
              "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status"]

SERVER_GRAPHS = {
    "tods1": [
        ("com-dblp",     [3, 5, 8, 10, 15]),
        ("com-youtube",  [3, 5, 8, 12, 16]),
    ],
    "tods2": [
        # collaboration / e-commerce
        ("com-amazon.ungraph",       [3, 5, 8]),
        ("com-dblp",                 [3, 5, 8, 10, 15]),
        # social — small to medium
        ("com-youtube",              [3, 5, 8, 12, 16]),
        ("twitter_combined",         [3, 5, 8, 12]),
        ("wiki-Talk",                [3, 5, 8, 10]),
        ("soc-pokec-relationships",  [3, 5, 8]),
        # com-lj (35M edges) and com-orkut (117M edges) were probed but their
        # SDCT leaf set exceeds the host's memory budget on a shared server:
        # com-lj s=3 → 335 GB peak RSS, com-lj s=5 → 305 GB peak, both
        # OOM-killed (probe_logs/com-lj_s{3,5}_T24.log). Excluded from the
        # sweep — the "very large graph" point in this experiment is
        # web-it-2004 (1.15B edges, sparse) which is tractable.
        # web
        ("web-Google",               [3, 5, 8, 12, 20]),
        ("web-Stanford",             [3, 5, 10, 20, 40, 60]),
        ("web-it-2004",              [3, 10, 30, 100, 200, 400]),
        # tiny sanity (kept as "no parallel benefit" reference cell)
        ("dblp-core30",              [3, 5, 10, 20]),
    ],
}

THREADS = [1, 2, 4, 8, 16, 24, 32, 48, 64]
RUNS_PER_CFG = 3   # take median across runs


_NUM = r"-?[\d.]+(?:e[+-]?\d+)?"
_PAR4_RE = re.compile(
    rf"PAR4 graph=(\S+) s=(\d+) T=(\d+) ms=(\d+) leaves=({_NUM})(?: ref_leaves=({_NUM}))? verify=(\S+)"
)


def parse_par4(txt: str) -> dict | None:
    for line in txt.splitlines():
        m = _PAR4_RE.search(line)
        if m:
            return {
                "graph": m.group(1), "s": int(m.group(2)),
                "threads": int(m.group(3)), "ms": int(m.group(4)),
                "leaves": float(m.group(5)),
                "ref_leaves": float(m.group(6)) if m.group(6) else -1.0,
                "verify": m.group(7),
            }
    return None


def load_existing_keys(csv_path: Path) -> set:
    """Count rows per (graph, s, threads); skip if we have ≥ RUNS_PER_CFG."""
    import csv as _csv
    counts = {}
    if not csv_path.exists(): return set()
    with csv_path.open() as f:
        for row in _csv.DictReader(f):
            if row.get("status") != "OK": continue
            try:
                k = (row["graph"], int(row["s"]), int(row["threads"]))
            except (ValueError, KeyError):
                continue
            counts[k] = counts.get(k, 0) + 1
    return {k for k, n in counts.items() if n >= RUNS_PER_CFG}


def append_row(csv_path: Path, row: dict) -> None:
    import csv as _csv
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("a", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=FIELDNAMES)
        if write_header: w.writeheader()
        w.writerow({k: row.get(k, "") for k in FIELDNAMES})


_TIME_BIN = "/usr/bin/time"


def parse_time_v(stderr: str) -> dict:
    """Parse the verbose `/usr/bin/time -v` block out of the binary's stderr.

    Recognises the labels emitted by GNU time-1.7+. Anything missing yields
    an empty string in the returned dict so the CSV column stays well-formed.
    """
    out = {k: "" for k in (
        "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
        "time_pagefaults_major", "time_pagefaults_minor",
        "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status")}
    for line in stderr.splitlines():
        line = line.strip()
        if line.startswith("Maximum resident set size"):
            out["time_max_rss_kB"]        = line.split(":")[-1].strip()
        elif line.startswith("User time (seconds)"):
            out["time_user_sec"]          = line.split(":")[-1].strip()
        elif line.startswith("System time (seconds)"):
            out["time_sys_sec"]           = line.split(":")[-1].strip()
        elif line.startswith("Elapsed (wall clock) time"):
            out["time_elapsed"]           = line.split(": ", 1)[-1].strip()
        elif line.startswith("Major (requiring I/O) page faults"):
            out["time_pagefaults_major"]  = line.split(":")[-1].strip()
        elif line.startswith("Minor (reclaiming a frame) page faults"):
            out["time_pagefaults_minor"]  = line.split(":")[-1].strip()
        elif line.startswith("Voluntary context switches"):
            out["time_voluntary_ctxt"]    = line.split(":")[-1].strip()
        elif line.startswith("Involuntary context switches"):
            out["time_involuntary_ctxt"]  = line.split(":")[-1].strip()
        elif line.startswith("Exit status"):
            out["time_exit_status"]       = line.split(":")[-1].strip()
    return out


def run_one(graph: str, s: int, T: int, run_idx: int, verify: bool) -> dict:
    gpath = f"graphs/{graph}.edges"
    log_path = LOGDIR / f"{graph}_s{s}_T{T}_r{run_idx}.log"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(T)
    # Pin threads to physical cores. tods2 is 2 sockets × 24 cores × 2 HT.
    # close+cores gives stable scaling up to 24 (one socket); past that NUMA
    # crossing is honest and reproducible.
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"]    = "cores"
    # /usr/bin/time -v wrapper: captures peak RSS, page faults, ctxt switches,
    # and exit status (incl. signals like SIGKILL=137 for OOM).  The verbose
    # report is written to stderr — we always persist combined stdout+stderr.
    args = [_TIME_BIN, "-v", BIN, gpath, str(s), str(T)]
    if verify: args.append("verify")
    cmd_line = " ".join(args)
    t0 = time.time()
    try:
        proc = subprocess.run(args, env=env, capture_output=True, text=True,
                              timeout=TIMEOUT)
        wall_ms = (time.time() - t0) * 1000.0
        # Always persist the full log, regardless of exit status.
        log_path.write_text(
            f"[CMD] {cmd_line}\n"
            f"[ENV] OMP_NUM_THREADS={T} OMP_PROC_BIND=close OMP_PLACES=cores\n"
            f"[RC]  {proc.returncode}\n\n"
            f"--STDOUT--\n{proc.stdout}\n"
            f"--STDERR--\n{proc.stderr}\n"
        )
        time_v = parse_time_v(proc.stderr)
        if proc.returncode == 0:
            par = parse_par4(proc.stdout)
            if par is None:
                return {"status": "PARSE_FAIL", "wall_ms": wall_ms, **time_v}
            return {"status": "OK", "wall_ms": wall_ms, "build_ms": par["ms"],
                    "leaves": par["leaves"], "ref_leaves": par["ref_leaves"],
                    "verify": par["verify"], **time_v}
        else:
            return {"status": f"ERROR({proc.returncode})",
                    "wall_ms": wall_ms, **time_v}
    except subprocess.TimeoutExpired as e:
        # Even on timeout, capture whatever output was produced so we can
        # see how far the run got and what its memory profile looked like.
        out = e.stdout.decode("utf-8", "replace") if e.stdout else ""
        err = e.stderr.decode("utf-8", "replace") if e.stderr else ""
        log_path.write_text(
            f"[CMD] {cmd_line}\n"
            f"[ENV] OMP_NUM_THREADS={T} OMP_PROC_BIND=close OMP_PLACES=cores\n"
            f"[TIMEOUT after {TIMEOUT}s]\n\n"
            f"--STDOUT--\n{out}\n"
            f"--STDERR--\n{err}\n"
        )
        return {"status": "TIMEOUT", "wall_ms": TIMEOUT * 1000.0,
                **parse_time_v(err)}


def main():
    print("=" * 60)
    print(f"  bench_par_sdct  {time.strftime('%F %T')}")
    print("=" * 60)
    cfg = ServerConfig.detect(DEFAULT_SERVERS)
    print(f"server: {cfg.name}", flush=True)

    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    build_main(["bench_sdct_par_scale"])

    graphs = SERVER_GRAPHS.get(cfg.name, SERVER_GRAPHS["tods2"])
    avail = link_graphs([g for g, _ in graphs], cfg)
    graphs = [(g, ss) for g, ss in graphs if g in avail]
    print(f"graphs: {[g for g, _ in graphs]}", flush=True)

    skip_full = load_existing_keys(OUTCSV)
    print(f"skipping {len(skip_full)} fully-done (graph,s,T) cells", flush=True)

    total = sum(len(ss) for _, ss in graphs) * len(THREADS) * RUNS_PER_CFG
    done_cnt = 0
    for graph, s_list in graphs:
        for s in s_list:
            # verify on T=1, run 0 only — others trust verify=OK from T=1
            for T in THREADS:
                if (graph, s, T) in skip_full:
                    done_cnt += RUNS_PER_CFG
                    continue
                for run_idx in range(RUNS_PER_CFG):
                    verify = (T == 1 and run_idx == 0)
                    r = run_one(graph, s, T, run_idx, verify)
                    row = {
                        "graph": graph, "s": s, "threads": T,
                        "status": r["status"],
                        "wall_ms":  f"{r['wall_ms']:.1f}",
                        "build_ms": r.get("build_ms", ""),
                        "leaves":   r.get("leaves", ""),
                        "ref_leaves": r.get("ref_leaves", ""),
                        "verify":   r.get("verify", ""),
                    }
                    append_row(OUTCSV, row)
                    done_cnt += 1
                    print(f"  [{done_cnt}/{total}] {graph} s={s} T={T} run={run_idx} "
                          f"{r['status']} build={r.get('build_ms', 'N/A')}ms",
                          flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
