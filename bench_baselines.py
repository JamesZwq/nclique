#!/usr/bin/env python3
"""
Experiment 6 — Comparison with non-CPI baselines.

The paper currently compares only against the CPI-based SOTA of
sigmodCPI (Reference). For a stronger evaluation we want at least one
non-CPI baseline:

    (a) Sariyuce et al., TODS 2022   — original (r,s)-nucleus decomposition
    (b) Fang et al., PVLDB 2019      — (k, Psi)-core (Psi=K_s)

Neither is currently in this codebase. The pinned tods2 home dir contains
~/decomp_baselines.cpp (per project memory) and ~/maxclq /
~/arb-nucleus-* — those are likely candidate sources.

This script builds a comparison harness around the NEW main binary
(REF + ST_V3) plus any external baseline binaries you supply via
BASELINE_BINS. Each baseline is invoked with (graph, r, s) as positional
args (extend `run_baseline` if a baseline takes a different CLI).

To activate a baseline, set BASELINE_BINS below and ensure the binary is
built and on $PATH.

Usage:
    nohup python3 bench_baselines.py tods2 > bench_baselines.log 2>&1 &
"""
from __future__ import annotations
import os, re, subprocess, time
from pathlib import Path

from bench_lib import (
    DEFAULT_SERVERS, ServerConfig, ParallelRunner, Job,
    raise_stack, link_graphs, load_done, build_main,
)

raise_stack()

BIN = "./build/bin/degeneracy_cliques"
TIMEOUT = 3600
OUTCSV = Path("paper_data/bench_baselines.csv")
LOGDIR = Path("bench_baselines_logs")
FIELDNAMES = ["graph", "r", "s", "algo", "status",
              "wall_ms", "total_ms", "peel_ms", "mem_kB"]

# (graph, [(r, s), ...]) — kept small until baselines are wired up.
SERVER_GRAPHS = {
    "tods1": [
        ("com-dblp",     [(1, 5), (1, 10), (1, 15)]),
        ("com-youtube",  [(1, 5), (1, 8), (1, 12)]),
    ],
    "tods2": [
        ("web-Stanford", [(1, 5), (1, 10), (1, 20), (1, 40)]),
        ("dblp-core30",  [(1, 5), (1, 10), (1, 20)]),
    ],
}

# ---- our entries (always available) ----
OUR_ALGOS = [
    ("REF",   {}),                                   # mutable-tree baseline (a.k.a. SOTA)
    ("ST_V3", {"PIVOTER_RUN_ST_V3": "1"}),
]

# ---- external baselines: set to None or path to binary ----
# Each binary should accept: <graph> <r> <s> and emit `time: <ms>` somewhere.
# Replace with actual paths once the baseline binaries are built.
BASELINE_BINS: dict[str, str | None] = {
    # "Sariyuce_TODS22": "~/maxclq/sariyuce_nucleus",
    # "Fang_PVLDB19":    "~/arb-nucleus-decomp/kpsi_core",
    "Sariyuce_TODS22": None,
    "Fang_PVLDB19":    None,
}


def parse_baseline_time(txt: str) -> float:
    """Best-effort regex; baselines vary in their timing format."""
    for pat in (r'time:\s*([\d.]+)\s*ms',
                r'total_time_ms\s*[:=]\s*([\d.]+)',
                r'elapsed:\s*([\d.]+)'):
        m = re.search(pat, txt, re.I)
        if m:
            try: return float(m.group(1))
            except: pass
    return -1.0


def gen_our_jobs(graphs, done):
    for graph, rs_list in graphs:
        gpath = f"graphs/{graph}.edges"
        if not os.path.exists(gpath): continue
        for r, s in rs_list:
            for algo, env in OUR_ALGOS:
                key = (graph, str(r), str(s), algo)
                if key in done: continue
                yield Job(
                    key=key,
                    cmd=[BIN, gpath, str(r), str(s)],
                    env=dict(env, OMP_NUM_THREADS="1"),
                    log_path=LOGDIR / f"{graph}_r{r}_s{s}_{algo}.log",
                    timeout=TIMEOUT,
                    extra={"graph": graph, "r": r, "s": s, "algo": algo,
                           "kind": "ours"},
                )


def gen_baseline_jobs(graphs, done):
    for algo, bin_path in BASELINE_BINS.items():
        if bin_path is None: continue
        bin_path = os.path.expanduser(bin_path)
        if not os.path.isfile(bin_path) or not os.access(bin_path, os.X_OK):
            print(f"  [skip] baseline {algo}: {bin_path} not executable", flush=True)
            continue
        for graph, rs_list in graphs:
            gpath = f"graphs/{graph}.edges"
            if not os.path.exists(gpath): continue
            for r, s in rs_list:
                key = (graph, str(r), str(s), algo)
                if key in done: continue
                yield Job(
                    key=key,
                    cmd=[bin_path, gpath, str(r), str(s)],
                    env={"OMP_NUM_THREADS": "1"},
                    log_path=LOGDIR / f"{graph}_r{r}_s{s}_{algo}.log",
                    timeout=TIMEOUT,
                    extra={"graph": graph, "r": r, "s": s, "algo": algo,
                           "kind": "baseline"},
                )


def main():
    print("=" * 60)
    print(f"  bench_baselines  {time.strftime('%F %T')}")
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
        if ex["kind"] == "baseline":
            t_ms = parse_baseline_time(log_text)
            row_total = f"{t_ms:.1f}" if t_ms >= 0 else ""
            row_peel = ""
        else:
            row_total = f"{parsed['total_ms']:.1f}" if parsed["total_ms"] >= 0 else ""
            row_peel  = f"{parsed['peel_ms']:.1f}"  if parsed["peel_ms"]  >= 0 else ""
        runner.append_row({
            "graph": ex["graph"], "r": ex["r"], "s": ex["s"], "algo": ex["algo"],
            "status": status,
            "wall_ms":  f"{parsed['wall_ms']:.1f}",
            "total_ms": row_total,
            "peel_ms":  row_peel,
            "mem_kB":   f"{parsed['mem_kB']:.0f}" if parsed["mem_kB"] >= 0 else "",
        })
        tag = f"{ex['graph']:>14} r={ex['r']} s={ex['s']:>3} {ex['algo']:>14}"
        print(f"  {tag} {status} wall={parsed['wall_ms']:.0f}ms total={row_total or 'N/A'}",
              flush=True)

    # Run our jobs first (always available), then baselines.
    print("\n=== Phase 1: our algos ===", flush=True)
    runner.run(gen_our_jobs(graphs, done), on_finish=on_finish)

    if any(b is not None for b in BASELINE_BINS.values()):
        print("\n=== Phase 2: external baselines ===", flush=True)
        done = load_done(OUTCSV, ("graph", "r", "s", "algo"))
        runner.run(gen_baseline_jobs(graphs, done), on_finish=on_finish)
    else:
        print("\n[skip] No external baselines configured. Edit BASELINE_BINS.",
              flush=True)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
