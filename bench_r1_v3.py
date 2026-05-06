#!/usr/bin/env python3
"""
r=1 main benchmark, V3 / Pure SOTA against REF, 3 runs per cell.

Replaces bench_r1_main.py for paper §7.2-§7.3 (Figure 5 endtoend +
ten-graph time/RSS panels).  The previous CSV (`01_main_benchmark_762`)
was produced with `PIVOTER_RUN_ST=1`, the FIRST static variant — paper
now claims "Pure" which means V3 (`PIVOTER_RUN_ST_V3=1`, with fused
SDCT + event-driven PeelH).  Same labelling discipline as bench_par_sdct:
each run wrapped in /usr/bin/time -v, full stdout+stderr persisted.

Outputs:
    paper_data/01_main_benchmark_v3.csv      # new CSV, 3 runs/cell
    bench_r1_v3_logs/<graph>_s<s>_<algo>_r<idx>.log

Usage:
    nohup python3 bench_r1_v3.py > /tmp/bench_r1_v3.log 2>&1 &
"""
from __future__ import annotations
import csv as _csv
import os, re, subprocess, sys, threading, time
from pathlib import Path

# Stack-limit raise (deep BK recursion at large s).
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
    print(f"[stack] soft={resource.getrlimit(resource.RLIMIT_STACK)[0]/1024/1024:.0f}MB",
          flush=True)
except Exception as e:
    print(f"[stack] {e}", flush=True)

BIN          = "./build/bin/degeneracy_cliques"
TIME_BIN     = "/usr/bin/time"
TIMEOUT      = 3600                          # per-cell ceiling (1 h)
RUNS_PER_CFG = 3
OUTCSV       = Path("paper_data/01_main_benchmark_v3.csv")
LOGDIR       = Path("bench_r1_v3_logs")

# Parallelism (memory-gated).  Each worker is one (graph, s, algo, run)
# subprocess. tods2 has 503 GB / 96 logical cores; web-it-2004 REF tops
# ~600 MB peak, twitter REF at high s tops ~30 GB.  Cap at MAX_WORKERS,
# and refuse to launch a new worker while used_mem > MEM_LIMIT_GB.
MAX_WORKERS    = int(os.environ.get("BENCH_WORKERS", "16"))
MEM_LIMIT_GB   = float(os.environ.get("BENCH_MEM_LIMIT_GB", "300"))
LAUNCH_POLL_S  = 3.0
SETTLE_AFTER_LAUNCH_S = 0.3


def get_used_mem_gb() -> float:
    try:
        info = {}
        with open("/proc/meminfo") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    info[parts[0].rstrip(":")] = int(parts[1])
        return (info.get("MemTotal", 0) - info.get("MemAvailable", 0)) / 1024 / 1024
    except Exception:
        return 0.0

# Per-graph s_max — same as the prior ST run so cells are 1:1 comparable.
GRAPHS: list[tuple[str, int]] = [
    ("com-dblp",                 30),
    ("com-amazon.ungraph",       30),
    ("twitter_combined",         30),
    ("web-Stanford",             61),
    ("web-Google",               40),
    ("com-youtube",              17),
    ("web-it-2004",              430),
    ("wiki-Talk",                25),
    ("soc-pokec-relationships",  25),
    ("com-orkut",                15),
]

ALGOS = [
    ("Pure",   {"PIVOTER_RUN_ST_V3": "1"}),  # V3 = fused SDCT + event-driven PeelH
    ("REF_R1", {}),                           # default: NCliqueVertexCoreDecomposition
]

FIELDNAMES = [
    "graph", "r", "s", "algorithm", "run", "status",
    "wall_ms", "took_ms", "build_ms", "peel_ms", "memory_kB",
    "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
    "time_pagefaults_major", "time_pagefaults_minor",
    "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status",
]
# Why both build_ms and peel_ms separately: paper §7.2 quotes peel-only
# speedup (V3 vs REF differ in peel phase only — both share the same
# SDCT_Fused build). took_ms = build + peel for legacy compatibility but
# the actual algorithm-vs-algorithm comparison must use peel_ms ratio.

# parse_timing greps the binary's own clock (peel + build); wall_ms is the
# subprocess.run measurement, max_rss_kB comes from /usr/bin/time -v.
_TOOK_RE = re.compile(r"NucleusCoreDecomposition took:\s*([\d.]+)\s*ms")
_TOTAL_OURS_RE = re.compile(r"ST_V3 r=1 \(peel\) took:\s*([\d.]+)\s*ms")
_BUILD_RE = re.compile(r"(?:SDCT_Fused|SDCT\+callback) took:\s*([\d.]+)\s*ms")
_MEM_RE   = re.compile(r"Final Memory:\s*([\d.]+)")


def parse_runtime(stdout: str) -> tuple[float, float, float]:
    """Return (build_ms, peel_ms, mem_kB) from binary stdout. -1 = missing."""
    build = -1.0; peel = -1.0; mem = -1.0
    if (m := _BUILD_RE.search(stdout)): build = float(m.group(1))
    # V3 prints both "ST_V3 r=1 (peel)" and an outer "NucleusCoreDecomposition";
    # REF prints only "NucleusCoreDecomposition".  Peel = whichever applies.
    if (m := _TOTAL_OURS_RE.search(stdout)): peel = float(m.group(1))
    elif (m := _TOOK_RE.search(stdout)):     peel = float(m.group(1))
    if (m := _MEM_RE.search(stdout)):        mem  = float(m.group(1))
    return build, peel, mem


def parse_time_v(stderr: str) -> dict:
    out = {k: "" for k in (
        "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
        "time_pagefaults_major", "time_pagefaults_minor",
        "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status")}
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


def load_existing_counts(csv_path: Path) -> dict:
    """Count OK rows per (graph, s, algo) for resume support."""
    counts: dict = {}
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


def append_row(csv_path: Path, row: dict) -> None:
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("a", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=FIELDNAMES)
        if write_header: w.writeheader()
        w.writerow({k: row.get(k, "") for k in FIELDNAMES})


def run_one(graph: str, s: int, algo: str, env_extra: dict, run_idx: int) -> dict:
    gpath = f"./graphs/{graph}.edges"
    log_path = LOGDIR / f"{graph}_s{s}_{algo}_r{run_idx}.log"
    env = os.environ.copy(); env.update(env_extra)
    args = [TIME_BIN, "-v", BIN, gpath, "1", str(s)]
    cmd_line = " ".join(args)
    t0 = time.time()
    try:
        proc = subprocess.run(args, env=env, capture_output=True, text=True,
                              timeout=TIMEOUT)
        wall_ms = (time.time() - t0) * 1000.0
        log_path.write_text(
            f"[CMD] {cmd_line}\n"
            f"[ENV] {env_extra}\n"
            f"[RC]  {proc.returncode}\n\n"
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
            f"[CMD] {cmd_line}\n[ENV] {env_extra}\n"
            f"[TIMEOUT after {TIMEOUT}s]\n\n"
            f"--STDOUT--\n{out}\n--STDERR--\n{err}\n"
        )
        return {"status": "TIMEOUT", "wall_ms": TIMEOUT * 1000.0,
                **parse_time_v(err)}


def main():
    print("=" * 60)
    print(f"  bench_r1_v3  {time.strftime('%F %T')}")
    print("=" * 60)
    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    # Verify graphs exist before kicking off (avoids minutes of wasted scheduling).
    avail = []
    for g, smax in GRAPHS:
        if Path(f"./graphs/{g}.edges").exists():
            avail.append((g, smax))
        else:
            print(f"[skip] {g} not in graphs/ — omitting", flush=True)
    print(f"graphs: {[g for g, _ in avail]}", flush=True)

    counts = load_existing_counts(OUTCSV)
    n_have_full = sum(1 for n in counts.values() if n >= RUNS_PER_CFG)
    print(f"already have ≥{RUNS_PER_CFG} OK runs for {n_have_full} (graph,s,algo) cells",
          flush=True)

    # Skip floor: once a (graph, algo) hits TIMEOUT or signal-9 OOM at some s,
    # don't keep grinding higher s for that algo.  Mutated by completion
    # callbacks; checked when scheduling.
    skip_floor: dict = {}
    csv_lock = threading.Lock()

    # Build the full job list as (graph, s, algo, env_extra, run_idx).
    jobs: list[tuple] = []
    for graph, s_max in avail:
        for s in range(2, s_max + 1):
            for algo, env_extra in ALGOS:
                already = counts.get((graph, s, algo), 0)
                for run_idx in range(already, RUNS_PER_CFG):
                    jobs.append((graph, s, algo, env_extra, run_idx))
    print(f"[plan] {len(jobs)} jobs to schedule across ≤{MAX_WORKERS} workers "
          f"(mem gate {MEM_LIMIT_GB:.0f} GB)", flush=True)

    done_cnt = 0
    total = len(jobs)

    def _emit(row: dict) -> None:
        with csv_lock:
            append_row(OUTCSV, row)

    def _job_callback(job, result):
        nonlocal done_cnt
        graph, s, algo, _, run_idx = job
        row = {
            "graph": graph, "r": 1, "s": s, "algorithm": algo,
            "run": run_idx, "status": result["status"],
            "wall_ms": f"{result['wall_ms']:.1f}",
            "took_ms":  result.get("took_ms",  ""),
            "build_ms": result.get("build_ms", ""),
            "peel_ms":  result.get("peel_ms",  ""),
            "memory_kB": result.get("memory_kB", ""),
            **{k: result.get(k, "") for k in (
                "time_max_rss_kB","time_user_sec","time_sys_sec",
                "time_elapsed","time_pagefaults_major",
                "time_pagefaults_minor","time_voluntary_ctxt",
                "time_involuntary_ctxt","time_exit_status")},
        }
        _emit(row)
        done_cnt += 1
        print(f"  [{done_cnt}/{total}] {graph} s={s} {algo} run={run_idx} "
              f"{result['status']} wall={result['wall_ms']/1000:.2f}s "
              f"rss={result.get('time_max_rss_kB','?')}kB",
              flush=True)
        if result["status"] == "TIMEOUT" or result["status"].startswith("FAIL(-9"):
            with csv_lock:
                cur = skip_floor.get((graph, algo))
                if cur is None or s < cur:
                    skip_floor[(graph, algo)] = s
                    print(f"    -> skip-floor: {graph}/{algo} s>={s}", flush=True)

    # Schedule via ThreadPoolExecutor: each worker spends ~all its time
    # waiting on subprocess.run, so threads (with the GIL released across
    # subprocess calls) are sufficient — no need for full multiprocessing.
    from concurrent.futures import ThreadPoolExecutor, FIRST_COMPLETED, wait

    print(f"[setup] threadpool max_workers={MAX_WORKERS}", flush=True)
    pool = ThreadPoolExecutor(max_workers=MAX_WORKERS)
    inflight: dict = {}  # future -> job tuple
    job_iter = iter(jobs)
    pending_job = next(job_iter, None)
    print(f"[setup] first pending: {pending_job}", flush=True)

    iter_count = 0
    last_progress_t = time.time()
    while pending_job is not None or inflight:
        # Schedule new jobs while there's room and memory headroom.
        while pending_job is not None and len(inflight) < MAX_WORKERS:
            graph, s, algo, env_extra, run_idx = pending_job
            # Skip-floor (this run-and-everything-higher-s for this algo).
            sf = skip_floor.get((graph, algo))
            if sf is not None and s >= sf:
                _emit({"graph": graph, "r": 1, "s": s, "algorithm": algo,
                       "run": run_idx, "status": "SKIP_FLOOR"})
                done_cnt += 1
                pending_job = next(job_iter, None)
                continue
            mem = get_used_mem_gb()
            if mem > MEM_LIMIT_GB:
                # Memory pressure — wait for someone to finish before launching.
                break
            future = pool.submit(run_one, graph, s, algo, env_extra, run_idx)
            inflight[future] = pending_job
            print(f"  [submit] {graph} s={s} {algo} run={run_idx} "
                  f"(inflight={len(inflight)})", flush=True)
            time.sleep(SETTLE_AFTER_LAUNCH_S)
            pending_job = next(job_iter, None)
        # Wait for at least one job to finish (or poll if nothing in flight).
        if not inflight:
            time.sleep(LAUNCH_POLL_S)
            continue
        # Heartbeat every ~60s while waiting.
        if time.time() - last_progress_t > 60:
            print(f"  [heartbeat] inflight={len(inflight)} done={done_cnt}/{total} "
                  f"mem={get_used_mem_gb():.1f}GB", flush=True)
            last_progress_t = time.time()
        done_set, _ = wait(list(inflight.keys()), timeout=LAUNCH_POLL_S,
                           return_when=FIRST_COMPLETED)
        for fut in done_set:
            job = inflight.pop(fut)
            try:
                result = fut.result()
            except Exception as exc:
                result = {"status": f"EXCEPT({exc!r})", "wall_ms": 0.0}
            _job_callback(job, result)
            last_progress_t = time.time()
        iter_count += 1
    pool.shutdown(wait=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
