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
import os, re, subprocess, sys, time
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
    "wall_ms", "took_ms", "memory_kB",
    "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
    "time_pagefaults_major", "time_pagefaults_minor",
    "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status",
]

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
                    "took_ms": took if took >= 0 else "",
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
    # don't keep grinding higher s for that algo.
    skip_floor: dict = {}

    total = sum((smax - 1) for _, smax in avail) * len(ALGOS) * RUNS_PER_CFG
    done_cnt = 0
    for graph, s_max in avail:
        for s in range(2, s_max + 1):
            for algo, env_extra in ALGOS:
                already = counts.get((graph, s, algo), 0)
                if already >= RUNS_PER_CFG:
                    done_cnt += RUNS_PER_CFG
                    continue
                if (graph, algo) in skip_floor and s >= skip_floor[(graph, algo)]:
                    # Skip-floor: write a SKIP row per missing run for clarity.
                    for run_idx in range(already, RUNS_PER_CFG):
                        append_row(OUTCSV, {
                            "graph": graph, "r": 1, "s": s, "algorithm": algo,
                            "run": run_idx, "status": "SKIP_FLOOR"})
                        done_cnt += 1
                    continue
                for run_idx in range(already, RUNS_PER_CFG):
                    r = run_one(graph, s, algo, env_extra, run_idx)
                    append_row(OUTCSV, {
                        "graph": graph, "r": 1, "s": s, "algorithm": algo,
                        "run": run_idx,
                        "status":  r["status"],
                        "wall_ms": f"{r['wall_ms']:.1f}",
                        "took_ms": r.get("took_ms", ""),
                        "memory_kB": r.get("memory_kB", ""),
                        **{k: r.get(k, "") for k in (
                            "time_max_rss_kB","time_user_sec","time_sys_sec",
                            "time_elapsed","time_pagefaults_major",
                            "time_pagefaults_minor","time_voluntary_ctxt",
                            "time_involuntary_ctxt","time_exit_status")},
                    })
                    done_cnt += 1
                    print(f"  [{done_cnt}/{total}] {graph} s={s} {algo} run={run_idx} "
                          f"{r['status']} wall={r['wall_ms']/1000:.2f}s "
                          f"rss={r.get('time_max_rss_kB','?')}kB",
                          flush=True)
                    if r["status"] in ("TIMEOUT",) or r["status"].startswith("FAIL(-9"):
                        skip_floor[(graph, algo)] = s
                        print(f"    -> skip-floor: {graph}/{algo} s>={s}", flush=True)
                        break
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
