#!/usr/bin/env python3
"""
remeasure_mem.py — single-stream median peak-memory re-measurement.

Re-runs a fixed (graph, algo, r, s) grid sequentially (NO concurrency, so
no memory-pressure jitter) REPS times and records the median peak memory
(`[Memory-*] Final Memory`, which is VmHWM on Linux / ru_maxrss on macOS).

Used to clean run-to-run allocator noise out of the memory figure on the
small, fast graphs (e.g. dblp-core30) without re-running the whole sweep.

Usage:
  python3 scripts/remeasure_mem.py \
      --graph dblp-core30 --algos RegNDC,V3LM_HIER \
      --r 4,5,6,7 --reps 3 --out bench_remeasure_dblp.csv
"""
from __future__ import annotations
import argparse, csv, json, os, re, statistics, subprocess, sys, time
from pathlib import Path

BIN = Path("./build/bin/degeneracy_cliques")
ALGO_ENV = {
    "RegNDC":      "PIVOTER_RUN_REGION_V3LM",
    "V3LM_HIER":   "PIVOTER_RUN_REGION_V3LM_HIER",
    "V3LM_NOCPI":  "PIVOTER_RUN_REGION_V3LM_NOCPI",
}
MEM_RE  = re.compile(r"\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)\s*kB")
WALL_RE = re.compile(r"Total time:\s*(\d+)\s*ms")


def maxclique(graph: str) -> int:
    j = Path("bench_v3_max_cliques.json")
    if j.exists():
        d = json.loads(j.read_text())
        if graph in d:
            return int(d[graph])
    raise SystemExit(f"max clique for {graph} not in bench_v3_max_cliques.json")


def run_once(graph: str, algo: str, r: int, s: int, timeout: int):
    env = os.environ.copy()
    env[ALGO_ENV[algo]] = "1"
    try:
        p = subprocess.run([str(BIN), f"graphs/{graph}.edges", str(r), str(s), "degen"],
                           env=env, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, None
    txt = p.stdout + p.stderr
    m = MEM_RE.search(txt); w = WALL_RE.search(txt)
    return (float(m.group(1)) if m else None,
            float(w.group(1)) if w else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--algos", default="RegNDC,V3LM_HIER")
    ap.add_argument("--r", default="4,5,6,7")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    graph = args.graph
    algos = [a.strip() for a in args.algos.split(",") if a.strip()]
    rs = [int(x) for x in args.r.split(",")]
    mc = maxclique(graph)

    out = Path(args.out)
    done = set()
    if out.exists():
        with out.open() as f:
            for row in csv.DictReader(f):
                done.add((row["graph"], row["algo"], int(row["r"]), int(row["s"])))
    new = not out.exists()
    fout = out.open("a", newline="")
    w = csv.writer(fout)
    if new:
        w.writerow(["graph", "algo", "r", "s", "reps", "mem_kB", "wall_ms_median"])

    total = 0
    t0 = time.time()
    for r in rs:
        for s in range(r + 1, mc + 1):
            for algo in algos:
                if (graph, algo, r, s) in done:
                    continue
                mems, walls = [], []
                ok = True
                for _ in range(args.reps):
                    m, wl = run_once(graph, algo, r, s, args.timeout)
                    if m is None:
                        ok = False; break
                    mems.append(m); walls.append(wl or -1)
                if not ok:
                    continue   # skip timeouts; leave the CSV row to the sweep
                med_mem = statistics.median(mems)
                med_wall = statistics.median([x for x in walls if x >= 0] or [-1])
                w.writerow([graph, algo, r, s, args.reps, f"{med_mem:.0f}", f"{med_wall:.0f}"])
                fout.flush()
                total += 1
                if total % 50 == 0:
                    print(f"  {total} cells, {time.time()-t0:.0f}s elapsed, "
                          f"last {algo} r={r} s={s} mem={med_mem/1024:.1f}MB", flush=True)
    fout.close()
    print(f"[remeasure] {total} cells written to {out} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
