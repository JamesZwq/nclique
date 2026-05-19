#!/usr/bin/env python3
"""
Extension run for web-it-2004 only: SPIN-star + CND past s=56.

The main bench (bench_r1_v3.py) stopped at s=56 for web-it-2004; the
SPIN comparison swept to s=430. This script extends the SPIN-star/CND
coverage on web-it-2004 only, appending into the same CSV so that
the headline figure has all three curves over the same x-range.

To keep wall time tractable on the local Mac (single-process, no
parallel workers), we run 1 run/cell instead of 3 -- the median we
already have for s=2..56 is the headline number; the extension is
visual coverage past that point.

Usage:
    nohup python3 bench_wi_extend.py > /tmp/bench_wi_extend.log 2>&1 &
"""
from __future__ import annotations
import csv as _csv
import os, re, subprocess, sys, time
from pathlib import Path

# Stack limit raise.
try:
    import resource
    _BIG = 1 << 30
    _soft, _hard = resource.getrlimit(resource.RLIMIT_STACK)
    _t = _BIG if _hard == resource.RLIM_INFINITY else min(_hard, _BIG)
    if _t != _soft:
        try:
            resource.setrlimit(resource.RLIMIT_STACK, (_t, _hard))
        except (ValueError, OSError):
            resource.setrlimit(resource.RLIMIT_STACK, (max(_soft, _t-4096), _hard))
except Exception:
    pass

BIN          = "./build/bin/degeneracy_cliques"
TIME_BIN     = "/usr/bin/time"
TIMEOUT      = 3600
RUNS_PER_CFG = 1  # extension only -- main bench already has 3 runs/cell for s<=56
OUTCSV       = Path("paper_data/01_main_benchmark_v3.csv")
LOGDIR       = Path("bench_wi_extend_logs")

GRAPH = "web-it-2004"
S_START = 57
S_MAX   = 200  # cap for tractability; SPIN's coverage goes to 430

ALGOS = [
    ("Pure",   {"PIVOTER_RUN_ST_V3": "1"}),
    ("REF_R1", {}),
]

FIELDNAMES = [
    "graph", "r", "s", "algorithm", "run", "status",
    "wall_ms", "took_ms", "build_ms", "peel_ms", "memory_kB",
    "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
    "time_pagefaults_major", "time_pagefaults_minor",
    "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status",
]

_TOOK_RE = re.compile(r"NucleusCoreDecomposition took:\s*([\d.]+)\s*ms")
_TOTAL_OURS_RE = re.compile(r"ST_V3 r=1 \(peel\) took:\s*([\d.]+)\s*ms")
_BUILD_RE = re.compile(r"(?:SDCT_Fused|SDCT\+callback) took:\s*([\d.]+)\s*ms")
_MEM_RE   = re.compile(r"Final Memory:\s*([\d.]+)")


def parse_runtime(stdout: str):
    build = peel = mem = -1.0
    if (m := _BUILD_RE.search(stdout)): build = float(m.group(1))
    if (m := _TOTAL_OURS_RE.search(stdout)): peel = float(m.group(1))
    elif (m := _TOOK_RE.search(stdout)):     peel = float(m.group(1))
    if (m := _MEM_RE.search(stdout)):        mem  = float(m.group(1))
    return build, peel, mem


def parse_time_v(stderr: str):
    out = {k: "" for k in (
        "time_max_rss_kB","time_user_sec","time_sys_sec","time_elapsed",
        "time_pagefaults_major","time_pagefaults_minor",
        "time_voluntary_ctxt","time_involuntary_ctxt","time_exit_status")}
    for line in stderr.splitlines():
        line = line.strip()
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


def load_existing_counts(csv_path: Path):
    counts = {}
    if not csv_path.exists(): return counts
    with csv_path.open() as f:
        for r in _csv.DictReader(f):
            if r.get("status") != "OK": continue
            try:
                k = (r["graph"], int(r["s"]), r["algorithm"])
                counts[k] = counts.get(k, 0) + 1
            except (KeyError, ValueError):
                pass
    return counts


def append_row(csv_path: Path, row: dict):
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("a", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=FIELDNAMES)
        if write_header: w.writeheader()
        w.writerow({k: row.get(k, "") for k in FIELDNAMES})


def run_one(s: int, algo: str, env_extra: dict, run_idx: int) -> dict:
    gpath = f"./graphs/{GRAPH}.edges"
    log_path = LOGDIR / f"{GRAPH}_s{s}_{algo}_r{run_idx}.log"
    env = os.environ.copy(); env.update(env_extra)
    args = [TIME_BIN, "-v", BIN, gpath, "1", str(s)]
    t0 = time.time()
    try:
        proc = subprocess.run(args, env=env, capture_output=True, text=True,
                              timeout=TIMEOUT)
        wall_ms = (time.time() - t0) * 1000.0
        log_path.write_text(
            f"[CMD] {' '.join(args)}\n[ENV] {env_extra}\n[RC]  {proc.returncode}\n\n"
            f"--STDOUT--\n{proc.stdout}\n--STDERR--\n{proc.stderr}\n"
        )
        time_v = parse_time_v(proc.stderr)
        if proc.returncode == 0:
            build, peel, mem = parse_runtime(proc.stdout)
            took = (build + peel) if (build >= 0 and peel >= 0) else (peel if peel >= 0 else -1)
            return {"status": "OK", "wall_ms": wall_ms,
                    "took_ms":  took  if took  >= 0 else "",
                    "build_ms": build if build >= 0 else "",
                    "peel_ms":  peel  if peel  >= 0 else "",
                    "memory_kB": mem if mem >= 0 else "",
                    **time_v}
        else:
            return {"status": f"FAIL({proc.returncode})",
                    "wall_ms": wall_ms, **time_v}
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode("utf-8", "replace") if e.stdout else ""
        err = e.stderr.decode("utf-8", "replace") if e.stderr else ""
        log_path.write_text(
            f"[CMD] {' '.join(args)}\n[ENV] {env_extra}\n"
            f"[TIMEOUT after {TIMEOUT}s]\n\n"
            f"--STDOUT--\n{out}\n--STDERR--\n{err}\n"
        )
        return {"status": "TIMEOUT", "wall_ms": TIMEOUT * 1000.0,
                **parse_time_v(err)}


def main():
    print("=" * 60, flush=True)
    print(f"  bench_wi_extend  {time.strftime('%F %T')}", flush=True)
    print(f"  graph={GRAPH}  s={S_START}..{S_MAX}  runs/cell={RUNS_PER_CFG}", flush=True)
    print("=" * 60, flush=True)
    LOGDIR.mkdir(exist_ok=True)

    counts = load_existing_counts(OUTCSV)
    timeout_floor = {}
    done_cnt = total = 0
    for s in range(S_START, S_MAX + 1):
        for algo, env_extra in ALGOS:
            need = max(0, RUNS_PER_CFG - counts.get((GRAPH, s, algo), 0))
            for run_idx in range(counts.get((GRAPH, s, algo), 0),
                                 counts.get((GRAPH, s, algo), 0) + need):
                total += 1

    print(f"[plan] {total} runs to do", flush=True)
    cnt = 0
    for s in range(S_START, S_MAX + 1):
        for algo, env_extra in ALGOS:
            if timeout_floor.get(algo) and s >= timeout_floor[algo]:
                append_row(OUTCSV, {"graph": GRAPH, "r": 1, "s": s,
                                    "algorithm": algo, "run": 0,
                                    "status": "SKIP_FLOOR"})
                continue
            already = counts.get((GRAPH, s, algo), 0)
            need = max(0, RUNS_PER_CFG - already)
            for run_idx in range(already, already + need):
                cnt += 1
                t0 = time.time()
                result = run_one(s, algo, env_extra, run_idx)
                elapsed = time.time() - t0
                row = {"graph": GRAPH, "r": 1, "s": s, "algorithm": algo,
                       "run": run_idx, "status": result["status"],
                       "wall_ms": f"{result['wall_ms']:.1f}",
                       "took_ms":  result.get("took_ms", ""),
                       "build_ms": result.get("build_ms", ""),
                       "peel_ms":  result.get("peel_ms", ""),
                       "memory_kB": result.get("memory_kB", ""),
                       **{k: result.get(k, "") for k in (
                           "time_max_rss_kB","time_user_sec","time_sys_sec",
                           "time_elapsed","time_pagefaults_major",
                           "time_pagefaults_minor","time_voluntary_ctxt",
                           "time_involuntary_ctxt","time_exit_status")}}
                append_row(OUTCSV, row)
                print(f"[{cnt}/{total}] s={s} {algo} {result['status']} "
                      f"wall={elapsed:.1f}s peel={result.get('peel_ms','?')}ms",
                      flush=True)
                if result["status"] == "TIMEOUT":
                    timeout_floor[algo] = s
                    print(f"  -> skip-floor: {algo} s>={s}", flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
