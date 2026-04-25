#!/usr/bin/env python3
"""
Main r=1 benchmark across ALL available graphs.

Runs Ours_ST (PIVOTER_RUN_ST=1) and REF_R1 (no env flag) on every graph in
GRAPHS, sweeping s from 2 to that graph's per-graph S_MAX.  Output CSV is
identical schema to benchmark_all_results.csv:

    graph,r,s,algorithm,time_ms,memory_kB,status

Resume-friendly: if a (graph, s, algorithm) row with status=OK already
exists in the CSV, the run is skipped.  Skip-on-timeout: when one run
times out for a particular (graph, algorithm) at s=t, all subsequent
s' > t for that (graph, algorithm) are skipped.

Usage:
    python3 bench_r1_main.py \
        --bin ./build/bin/degeneracy_cliques \
        --graph-dir ./graphs \
        --out paper_data/01_main_benchmark_all_graphs.csv \
        --runs 1
"""

from __future__ import annotations

import argparse, csv, os, re, signal, subprocess, sys, time
from pathlib import Path

# Stack-limit raise (deep BK recursion at large s).
try:
    import resource
    _BIG = 1 << 30
    _soft, _hard = resource.getrlimit(resource.RLIMIT_STACK)
    if _hard == resource.RLIM_INFINITY:
        _t = _BIG
    elif _hard >= _BIG:
        _t = _BIG
    else:
        _t = _hard
    if _t != _soft:
        try:
            resource.setrlimit(resource.RLIMIT_STACK, (_t, _hard))
        except (ValueError, OSError):
            resource.setrlimit(resource.RLIMIT_STACK, (max(_soft, _t-4096), _hard))
    print(f"[stack] soft={resource.getrlimit(resource.RLIMIT_STACK)[0]/1024/1024:.0f}MB",
          flush=True)
except Exception as e:
    print(f"[stack] {e}", flush=True)


# ============ Per-graph s ranges ============
# (graph_stem, s_max).  s starts at 2 always.
GRAPHS = [
    ("com-dblp",                30),
    ("com-amazon.ungraph",      30),
    ("twitter_combined",        30),
    ("web-Stanford",            61),
    ("web-Google",              40),
    ("com-youtube",             17),
    ("web-it-2004",            430),
    ("wiki-Talk",               25),
    ("soc-pokec-relationships", 25),
    ("com-orkut",               15),
]

ALGOS = [
    ("Ours_ST", {"PIVOTER_RUN_ST": "1"}),
    ("REF_R1",  {}),
]


def parse_existing(csv_path: Path) -> set:
    """Return {(graph, s, algo)} of rows already recorded with status==OK."""
    done = set()
    if not csv_path.exists():
        return done
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            if int(row.get("r", "0")) != 1:
                continue
            if row.get("status") != "OK":
                continue
            done.add((row["graph"], int(row["s"]), row["algorithm"]))
    return done


def parse_timing(stdout: str) -> tuple[float, float]:
    """Extract (total_ms, mem_kb) from binary stdout."""
    m_total = re.search(r'NucleusCoreDecomposition took:\s*([\d.]+)', stdout)
    m_mem   = re.search(r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)', stdout)
    total = float(m_total.group(1)) if m_total else -1.0
    mem   = float(m_mem.group(1))   if m_mem   else -1.0
    return total, mem


def run_one(bin_path: Path, gpath: Path, s: int, algo: str, env_extra: dict,
            timeout: int) -> tuple[str, float, float]:
    env = os.environ.copy(); env.update(env_extra)
    try:
        proc = subprocess.run(
            [str(bin_path), str(gpath), "1", str(s)],
            env=env, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return ("TIMEOUT", -1.0, -1.0)
    if proc.returncode != 0:
        return (f"FAIL({proc.returncode})", -1.0, -1.0)
    total_ms, mem_kb = parse_timing(proc.stdout)
    if total_ms < 0:
        return ("PARSE_FAIL", -1.0, -1.0)
    return ("OK", total_ms, mem_kb)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin",       default="./build/bin/degeneracy_cliques")
    ap.add_argument("--graph-dir", default="./graphs")
    ap.add_argument("--out",       default="bench_r1_all.csv")
    ap.add_argument("--runs",      type=int, default=1,
                    help="how many times to run each (graph, s, algo) — extra runs append fresh rows")
    ap.add_argument("--timeout",   type=int, default=3600,
                    help="per-run timeout in seconds (default 3600 = 1 hour)")
    args = ap.parse_args()

    bin_path = Path(args.bin); graph_dir = Path(args.graph_dir)
    out_csv = Path(args.out); out_csv.parent.mkdir(parents=True, exist_ok=True)
    if not bin_path.exists(): sys.exit(f"binary not found: {bin_path}")

    done = parse_existing(out_csv)
    write_header = not out_csv.exists() or out_csv.stat().st_size == 0

    fout = out_csv.open("a")
    if write_header:
        fout.write("graph,r,s,algorithm,time_ms,memory_kB,status\n")
        fout.flush()

    n_total = 0; n_done = 0; n_skip = 0; n_fail = 0
    # Track timeout floor per (graph, algo) within this session.
    timeout_floor = {}

    for graph_stem, s_max in GRAPHS:
        gpath = graph_dir / f"{graph_stem}.edges"
        if not gpath.exists():
            print(f"[skip] {graph_stem}: file not found", flush=True)
            continue

        for s in range(2, s_max + 1):
            for algo_label, env_extra in ALGOS:
                key = (graph_stem, s, algo_label)
                key_done_count = sum(1 for (g, ss, a) in done
                                     if (g, ss, a) == key)
                # Use simple "first run already done" check; the script's
                # main resume target is partial benchmarks, not exact replay.
                if (graph_stem, s, algo_label) in done:
                    n_skip += 1
                    continue

                floor = timeout_floor.get((graph_stem, algo_label))
                if floor is not None and s >= floor:
                    fout.write(f"{graph_stem},1,{s},{algo_label},,,SKIP_TIMEOUT_FLOOR\n")
                    fout.flush()
                    n_skip += 1
                    continue

                n_total += 1
                t0 = time.time()
                status, total_ms, mem_kb = run_one(
                    bin_path, gpath, s, algo_label, env_extra, args.timeout)
                elapsed = time.time() - t0
                t_str = f"{total_ms:.1f}" if total_ms >= 0 else ""
                m_str = f"{mem_kb:.0f}"   if mem_kb   >= 0 else ""
                fout.write(f"{graph_stem},1,{s},{algo_label},{t_str},{m_str},{status}\n")
                fout.flush()

                tag = f"[{graph_stem} s={s} {algo_label}]"
                print(f"{tag} {status} t={total_ms:.0f}ms rss={mem_kb:.0f}kB "
                      f"(wall={elapsed:.1f}s)", flush=True)

                if status == "TIMEOUT":
                    timeout_floor[(graph_stem, algo_label)] = s
                    n_fail += 1
                    print(f"  -> skip-floor: {graph_stem}/{algo_label} s>={s}",
                          flush=True)
                elif status != "OK":
                    n_fail += 1
                else:
                    n_done += 1

    fout.close()
    print(f"\n=== summary ===")
    print(f"  done={n_done}  failed={n_fail}  skipped={n_skip}  total_attempted={n_total}")


if __name__ == "__main__":
    main()
