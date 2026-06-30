#!/usr/bin/env python3
"""
LazyPop vs V3 head-to-head r=1 (1,s)-nucleus benchmark.

Runs both algorithms serially with /usr/bin/time -v, captures wall time
breakdown (build, peel, total) and peak RSS. Three trials per cell, median.

Schema:
    graph,s,algo,trial,status,total_ms,build_ms,peel_ms,peak_rss_kB,
    pop_refresh,pop_bounce,note

Resume-friendly: skips cells with >= RUNS_PER_CFG OK rows already in CSV.

Usage on server:
    nohup python3 bench_lazypop.py tods2 > bench_lazypop.log 2>&1 &
"""
from __future__ import annotations
import csv, os, re, subprocess, sys, time
from pathlib import Path

from bench_lib import DEFAULT_SERVERS, ServerConfig, raise_stack, link_graphs, build_main

raise_stack()

BIN = "./build/bin/degeneracy_cliques"
TIMEOUT = 1800
OUTCSV = Path("paper_data/bench_lazypop_vs_v3.csv")
LOGDIR = Path("bench_lazypop_logs")
RUNS_PER_CFG = 3

FIELDS = ["graph", "s", "algo", "trial", "status",
          "total_ms", "build_ms", "peel_ms", "peel_us", "peak_rss_kB",
          "pop_refresh", "pop_bounce", "note"]

# (graph_stem, [s values to test]) per server.
# Pick graphs that exercise different leaf-count / max-core regimes so we
# can see when LazyPop wins big vs ties V3.
SERVER_GRAPHS = {
    "tods1": [
        ("com-dblp",                 [3, 5, 8, 10, 15, 20]),
        ("com-youtube",              [3, 5, 8, 12, 16]),
        ("com-orkut",                [3, 5, 7]),
        ("soc-pokec-relationships",  [3, 5, 8]),
        ("web-Stanford",             [3, 5, 10, 20, 40, 60]),
    ],
    "tods2": [
        ("com-amazon.ungraph",       [3, 5, 8]),
        ("com-dblp",                 [3, 5, 8, 10, 15, 20]),
        ("com-youtube",              [3, 5, 8, 12, 16]),
        ("twitter_combined",         [3, 5, 8, 12]),
        ("wiki-Talk",                [3, 5, 8, 10]),
        ("soc-pokec-relationships",  [3, 5, 8]),
        ("web-Google",               [3, 5, 8, 12, 20]),
        ("web-Stanford",             [3, 5, 10, 20, 40, 60]),
        ("web-it-2004",              [3, 10, 30, 100, 200]),
    ],
}

ALGOS = [
    ("V3",      "PIVOTER_RUN_ST_V3"),
    ("LazyPop", "PIVOTER_RUN_LAZYPOP"),
]


_NUM = r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'
_RE_TOTAL_NCD = re.compile(rf'NucleusCoreDecomposition took:\s*({_NUM})')
_RE_V3_BUILD  = re.compile(rf'ST_V3 Build took:\s*({_NUM})\s*ms')
_RE_V3_PEEL   = re.compile(rf'ST_V[23]\s+r=1\s+\(peel\)\s+took:\s*({_NUM})\s*ms')
_RE_V3_US     = re.compile(rf'STV3_PEEL_US:\s*(\d+)')
_RE_LZ_LINE   = re.compile(
    rf'LazyPop:\s+(?:Lean|V3)\s+Build\s+({_NUM})\s*ms,\s*peel\s+({_NUM})\s*ms\s*'
    rf'\(pop_refresh=(\d+),\s*pop_bounce=(\d+)\)'
)
_RE_LZ_TOTAL  = re.compile(rf'LazyPop r=1 took:\s*({_NUM})')
_RE_LZ_BUILD_US = re.compile(rf'LAZYPOP_BUILD_US:\s*(\d+)')
_RE_LZ_PEEL_US  = re.compile(rf'LAZYPOP_PEEL_US:\s*(\d+)')
_RE_RSS = re.compile(rf'Maximum resident set size[^:]*:\s*({_NUM})')

def parse(stdout: str, stderr: str, algo: str) -> dict:
    out = {"total_ms": "", "build_ms": "", "peel_ms": "", "peel_us": "",
           "peak_rss_kB": "", "pop_refresh": "", "pop_bounce": ""}
    txt = stdout + "\n" + stderr
    if algo == "V3":
        if (m := _RE_V3_BUILD.search(txt)): out["build_ms"] = m.group(1)
        if (m := _RE_V3_PEEL.search(txt)):  out["peel_ms"]  = m.group(1)
        if (m := _RE_V3_US.search(txt)):    out["peel_us"]  = m.group(1)
        if (m := _RE_TOTAL_NCD.search(txt)): out["total_ms"] = m.group(1)
    else:
        if (m := _RE_LZ_LINE.search(txt)):
            out["build_ms"]    = m.group(1)
            out["peel_ms"]     = m.group(2)
            out["pop_refresh"] = m.group(3)
            out["pop_bounce"]  = m.group(4)
        if (m := _RE_LZ_BUILD_US.search(txt)):
            # Override ms build with μs-derived (more precise) if present.
            out["build_ms"] = f"{int(m.group(1)) / 1000:.3f}"
        if (m := _RE_LZ_PEEL_US.search(txt)):
            out["peel_us"] = m.group(1)
        if (m := _RE_LZ_TOTAL.search(txt)): out["total_ms"] = m.group(1)
    if (m := _RE_RSS.search(stderr)): out["peak_rss_kB"] = m.group(1)
    return out

def load_done(csv_path: Path) -> set[tuple]:
    counts = {}
    if not csv_path.exists(): return set()
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            if row.get("status") != "OK": continue
            try:
                k = (row["graph"], int(row["s"]), row["algo"])
            except (KeyError, ValueError):
                continue
            counts[k] = counts.get(k, 0) + 1
    return {k for k, n in counts.items() if n >= RUNS_PER_CFG}

def append(csv_path: Path, row: dict) -> None:
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if write_header: w.writeheader()
        w.writerow({k: row.get(k, "") for k in FIELDS})

def run_one(graph: str, s: int, algo: str, env_flag: str, trial: int) -> dict:
    gpath = f"graphs/{graph}.edges"
    log = LOGDIR / f"{graph}_s{s}_{algo}_t{trial}.log"
    env = os.environ.copy()
    env[env_flag] = "1"
    env["OMP_NUM_THREADS"] = "1"
    cmd = ["/usr/bin/time", "-v", BIN, gpath, "1", str(s)]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=TIMEOUT)
        wall = (time.time() - t0) * 1000.0
        log.write_text(
            f"[CMD] {' '.join(cmd)}\n[ENV] {env_flag}=1 OMP_NUM_THREADS=1\n"
            f"[RC] {proc.returncode}\n\n--STDOUT--\n{proc.stdout}\n--STDERR--\n{proc.stderr}\n"
        )
        if proc.returncode != 0:
            return {"status": f"ERR({proc.returncode})", "note": "nonzero rc"}
        parsed = parse(proc.stdout, proc.stderr, algo)
        return {"status": "OK", **parsed}
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode("utf-8", "replace") if e.stdout else ""
        err = e.stderr.decode("utf-8", "replace") if e.stderr else ""
        log.write_text(f"[CMD] {' '.join(cmd)}\n[TIMEOUT {TIMEOUT}s]\n\n--STDOUT--\n{out}\n--STDERR--\n{err}\n")
        return {"status": "TIMEOUT", "note": f"timeout {TIMEOUT}s"}

def main():
    print("=" * 60)
    print(f"  bench_lazypop  {time.strftime('%F %T')}")
    print("=" * 60)
    cfg = ServerConfig.detect(DEFAULT_SERVERS)
    print(f"server: {cfg.name}", flush=True)

    LOGDIR.mkdir(exist_ok=True)
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)

    build_main(["degeneracy_cliques"])

    graphs = SERVER_GRAPHS.get(cfg.name, SERVER_GRAPHS["tods2"])
    avail = link_graphs([g for g, _ in graphs], cfg)
    graphs = [(g, ss) for g, ss in graphs if g in avail]
    print(f"graphs: {[g for g,_ in graphs]}", flush=True)

    done = load_done(OUTCSV)
    print(f"already done: {len(done)} (graph,s,algo) cells", flush=True)

    total = sum(len(ss) for _, ss in graphs) * len(ALGOS) * RUNS_PER_CFG
    n = 0
    for graph, s_list in graphs:
        for s in s_list:
            for algo, env_flag in ALGOS:
                if (graph, s, algo) in done:
                    n += RUNS_PER_CFG
                    continue
                for trial in range(RUNS_PER_CFG):
                    r = run_one(graph, s, algo, env_flag, trial)
                    append(OUTCSV, {
                        "graph": graph, "s": s, "algo": algo, "trial": trial,
                        **r
                    })
                    n += 1
                    print(f"  [{n}/{total}] {graph} s={s} {algo} t={trial}: "
                          f"{r.get('status')} total={r.get('total_ms','?')} "
                          f"build={r.get('build_ms','?')} peel={r.get('peel_ms','?')}ms"
                          f"/{r.get('peel_us','?')}us "
                          f"rss={r.get('peak_rss_kB','?')}",
                          flush=True)
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
