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
              "leaves", "ref_leaves", "verify"]

SERVER_GRAPHS = {
    "tods1": [
        ("com-dblp",     [3, 5, 8, 10, 15]),
        ("com-youtube",  [3, 5, 8, 12, 16]),
    ],
    "tods2": [
        ("web-Stanford", [3, 5, 10, 20, 40, 60]),
        ("dblp-core30",  [3, 5, 10, 20]),
        ("web-it-2004",  [3, 10, 30, 100, 200, 400]),
    ],
}

THREADS = [1, 2, 4, 8, 16, 24]
RUNS_PER_CFG = 3   # take median across runs


_PAR4_RE = re.compile(
    r"PAR4 graph=(\S+) s=(\d+) T=(\d+) ms=(\d+) leaves=(\d+)(?: ref_leaves=(-?\d+))? verify=(\S+)"
)


def parse_par4(txt: str) -> dict | None:
    for line in txt.splitlines():
        m = _PAR4_RE.search(line)
        if m:
            return {
                "graph": m.group(1), "s": int(m.group(2)),
                "threads": int(m.group(3)), "ms": int(m.group(4)),
                "leaves": int(m.group(5)),
                "ref_leaves": int(m.group(6)) if m.group(6) else -1,
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


def run_one(graph: str, s: int, T: int, run_idx: int, verify: bool) -> dict:
    gpath = f"graphs/{graph}.edges"
    log_path = LOGDIR / f"{graph}_s{s}_T{T}_r{run_idx}.log"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(T)
    args = [BIN, gpath, str(s), str(T)]
    if verify: args.append("verify")
    t0 = time.time()
    try:
        proc = subprocess.run(args, env=env, capture_output=True, text=True,
                              timeout=TIMEOUT)
        wall_ms = (time.time() - t0) * 1000.0
        log_path.write_text(proc.stdout + "\n--STDERR--\n" + proc.stderr)
        if proc.returncode == 0:
            par = parse_par4(proc.stdout)
            if par is None:
                return {"status": "PARSE_FAIL", "wall_ms": wall_ms}
            return {"status": "OK", "wall_ms": wall_ms, "build_ms": par["ms"],
                    "leaves": par["leaves"], "ref_leaves": par["ref_leaves"],
                    "verify": par["verify"]}
        else:
            return {"status": f"ERROR({proc.returncode})", "wall_ms": wall_ms}
    except subprocess.TimeoutExpired:
        return {"status": "TIMEOUT", "wall_ms": TIMEOUT * 1000.0}


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
