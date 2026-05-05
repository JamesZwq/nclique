#!/usr/bin/env python3
"""
Phase-breakdown benchmark for paper §7.4 (Table tab:bd-time / tab:bd-mem,
Figure fig_phase_breakdown), V3 / Pure SOTA against REF.

Unlike bench_r1_v3.py (which aggregates total wall + total RSS for the
end-to-end speedup figure), this bench captures individual phase times
— load, build, peel — that the binary already prints, so the paper can
attribute the speedup to specific phases.

Replaces `02_breakdown_summary.csv` whose `algo="ours"` rows were
actually `PIVOTER_RUN_ST=1` (the older ST variant).  Output is a new
CSV; the legacy file stays for cross-checking.

Outputs:
    paper_data/02_breakdown_summary_v3.csv
    bench_r1_breakdown_v3_logs/<graph>_s<s>_<algo>_r<idx>.log

Usage:
    nohup python3 bench_r1_breakdown_v3.py > /tmp/bench_breakdown_v3.log 2>&1 &
"""
from __future__ import annotations
import csv as _csv
import os, re, subprocess, sys, time
from pathlib import Path

# Stack-limit raise.
try:
    import resource
    _BIG = 1 << 30
    _soft, _hard = resource.getrlimit(resource.RLIMIT_STACK)
    _t = _BIG if _hard == resource.RLIM_INFINITY else min(_hard, _BIG)
    if _t != _soft:
        try: resource.setrlimit(resource.RLIMIT_STACK, (_t, _hard))
        except (ValueError, OSError):
            resource.setrlimit(resource.RLIMIT_STACK, (max(_soft, _t-4096), _hard))
except Exception: pass

BIN          = "./build/bin/degeneracy_cliques"
TIME_BIN     = "/usr/bin/time"
TIMEOUT      = 3600
RUNS_PER_CFG = 3
OUTCSV       = Path("paper_data/02_breakdown_summary_v3.csv")
LOGDIR       = Path("bench_r1_breakdown_v3_logs")

# Same coverage as the legacy 02_breakdown_summary.csv so cells are 1:1
# diffable with the previous ST-based experiment.
GRAPHS_S = [
    ("com-youtube",    [3, 5, 8, 12, 16]),
    ("web-Stanford",   [3, 5, 10, 20, 40, 60]),
    ("web-it-2004",    [3, 10, 30, 100, 200, 400]),
]

ALGOS = [
    ("Pure",   {"PIVOTER_RUN_ST_V3": "1"}),
    ("REF_R1", {}),
]

FIELDNAMES = [
    "graph", "s", "algorithm", "run", "status",
    "wall_ms", "load_ms", "build_ms", "peel_ms",
    "total_time_ms", "memory_kB",
    "time_max_rss_kB", "time_user_sec", "time_sys_sec", "time_elapsed",
    "time_pagefaults_major", "time_pagefaults_minor",
    "time_voluntary_ctxt", "time_involuntary_ctxt", "time_exit_status",
]

# Phase regex.  Names match what daf::timeCount prints from the binary.
_LOAD_RE   = re.compile(r"loadAndSort took:\s*([\d.]+)\s*ms")
_BUILD_RE  = re.compile(r"(?:buildSDCT|SDCT_Fused|SDCT\+callback) took:\s*([\d.]+)\s*ms")
_PEEL_V3_RE = re.compile(r"ST_V3 r=1 \(peel\) took:\s*([\d.]+)\s*ms")
_PEEL_REF_RE = re.compile(r"NucleusCoreDecomposition took:\s*([\d.]+)\s*ms")
_MEM_RE    = re.compile(r"Final Memory:\s*([\d.]+)")


def parse_phases(stdout: str, algo: str) -> dict:
    """Extract per-phase times from binary stdout."""
    out = {"load_ms": "", "build_ms": "", "peel_ms": "", "memory_kB": ""}
    if (m := _LOAD_RE.search(stdout)):  out["load_ms"]  = float(m.group(1))
    if (m := _BUILD_RE.search(stdout)): out["build_ms"] = float(m.group(1))
    # V3 prints both ST_V3 (inner) and NucleusCoreDecomposition (outer).
    # REF prints only NucleusCoreDecomposition. Pick the inner V3 one when
    # available (more accurate peel timing).
    if algo == "Pure":
        if (m := _PEEL_V3_RE.search(stdout)):  out["peel_ms"] = float(m.group(1))
        elif (m := _PEEL_REF_RE.search(stdout)): out["peel_ms"] = float(m.group(1))
    else:
        if (m := _PEEL_REF_RE.search(stdout)): out["peel_ms"] = float(m.group(1))
    if (m := _MEM_RE.search(stdout)):  out["memory_kB"] = float(m.group(1))
    return out


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


def load_existing_counts(csv_path: Path) -> dict:
    counts: dict = {}
    if not csv_path.exists(): return counts
    with csv_path.open() as f:
        for r in _csv.DictReader(f):
            if r.get("status") != "OK": continue
            try:
                k = (r["graph"], int(r["s"]), r["algorithm"])
                counts[k] = counts.get(k, 0) + 1
            except (KeyError, ValueError): pass
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
            f"[CMD] {cmd_line}\n[ENV] {env_extra}\n"
            f"[RC]  {proc.returncode}\n\n"
            f"--STDOUT--\n{proc.stdout}\n--STDERR--\n{proc.stderr}\n"
        )
        time_v = parse_time_v(proc.stderr)
        if proc.returncode == 0:
            ph = parse_phases(proc.stdout, algo)
            total = ""
            try:
                total = (float(ph.get("load_ms") or 0) +
                         float(ph.get("build_ms") or 0) +
                         float(ph.get("peel_ms") or 0))
            except Exception:
                total = ""
            return {"status": "OK", "wall_ms": wall_ms,
                    **ph, "total_time_ms": total, **time_v}
        else:
            return {"status": f"FAIL({proc.returncode})",
                    "wall_ms": wall_ms, **time_v}
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode("utf-8","replace") if e.stdout else ""
        err = e.stderr.decode("utf-8","replace") if e.stderr else ""
        log_path.write_text(
            f"[CMD] {cmd_line}\n[ENV] {env_extra}\n"
            f"[TIMEOUT after {TIMEOUT}s]\n\n"
            f"--STDOUT--\n{out}\n--STDERR--\n{err}\n"
        )
        return {"status": "TIMEOUT", "wall_ms": TIMEOUT * 1000.0,
                **parse_time_v(err)}


def main():
    print("=" * 60)
    print(f"  bench_r1_breakdown_v3  {time.strftime('%F %T')}")
    print("=" * 60)
    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    avail = []
    for g, s_list in GRAPHS_S:
        if Path(f"./graphs/{g}.edges").exists():
            avail.append((g, s_list))
        else:
            print(f"[skip] {g} not in graphs/", flush=True)
    counts = load_existing_counts(OUTCSV)
    n_have_full = sum(1 for n in counts.values() if n >= RUNS_PER_CFG)
    print(f"already have ≥{RUNS_PER_CFG} OK runs for {n_have_full} (graph,s,algo) cells",
          flush=True)

    total = sum(len(ss) for _, ss in avail) * len(ALGOS) * RUNS_PER_CFG
    done_cnt = 0
    for graph, s_list in avail:
        for s in s_list:
            for algo, env_extra in ALGOS:
                already = counts.get((graph, s, algo), 0)
                if already >= RUNS_PER_CFG:
                    done_cnt += RUNS_PER_CFG
                    continue
                for run_idx in range(already, RUNS_PER_CFG):
                    r = run_one(graph, s, algo, env_extra, run_idx)
                    append_row(OUTCSV, {
                        "graph": graph, "s": s, "algorithm": algo,
                        "run": run_idx, "status": r["status"],
                        "wall_ms": f"{r['wall_ms']:.1f}",
                        "load_ms":   r.get("load_ms", ""),
                        "build_ms":  r.get("build_ms", ""),
                        "peel_ms":   r.get("peel_ms", ""),
                        "total_time_ms": r.get("total_time_ms", ""),
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
                          f"load={r.get('load_ms','?')} build={r.get('build_ms','?')} "
                          f"peel={r.get('peel_ms','?')}",
                          flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
