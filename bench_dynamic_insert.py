#!/usr/bin/env python3
"""
Correctness + performance driver for the dynamic (1,s)-core insert prototype.

Per sampled edge e=(u,v):
  coreBase = V3(G-e)          (the maintained state before the insert)
  run  dynamic_1s_core G-e s coreBase.tsv u v   -> CHANGED rows + STATS
  merged = coreBase patched with CHANGED
  VERIFY: merged == V3(G) per-vertex (missing id = core 0)

Baseline for speed: V3 static recompute on G, build+peel time only
(parsed from 'ST_V3 Build took' + 'STV3_PEEL_US', excludes load/sort —
conservative in the baseline's favor).

Usage:
  python3 bench_dynamic_insert.py --graph graphs/com-dblp.edges --s 5 \
      --edges 300 --workers 6 --tmpdir <scratch>
"""
from __future__ import annotations
import argparse, csv, os, random, re, subprocess, tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from bench_dynamic_locality import read_edge_file, run_v3

PROTO = os.environ.get("DYN1S_PROTO", "./build/bin/dynamic_1s_core")

_W = {}


def _init_worker(graph_path, s, workdir, full_cores):
    n_vertices, data_lines = read_edge_file(graph_path)
    _W.update(n=n_vertices, lines=data_lines, s=s,
              workdir=workdir, full=full_cores)


def one_edge(drop_idx: int):
    lines, n, s, workdir = _W["lines"], _W["n"], _W["s"], _W["workdir"]
    parts = lines[drop_idx].split()
    u, v = int(parts[0]), int(parts[1])
    tag = f"e{drop_idx}_{os.getpid()}"
    gfile = os.path.join(workdir, f"g_{tag}.edges")
    cfile = os.path.join(workdir, f"cb_{tag}.tsv")
    with open(gfile, "w") as f:
        f.write(f"{n} {len(lines) - 1}\n")
        for i, line in enumerate(lines):
            if i != drop_idx:
                f.write(line + "\n")
    try:
        base = run_v3(gfile, s, tag, workdir)
        with open(cfile, "w") as f:
            for x, c in base.items():
                f.write(f"{x}\t{c:.0f}\n")
        try:
            proc = subprocess.run([PROTO, gfile, str(s), cfile, str(u), str(v)],
                                  capture_output=True, text=True, timeout=600)
        except subprocess.TimeoutExpired:
            return dict(u=u, v=v, error="prototype timeout 600s")
        if proc.returncode != 0:
            return dict(u=u, v=v, error=proc.stderr[-200:])
    finally:
        os.unlink(gfile)
        if os.path.exists(cfile):
            os.unlink(cfile)

    merged = dict(base)
    stats = {}
    fallback = 0
    for line in proc.stdout.splitlines():
        if line.startswith("CHANGED"):
            _, x, old, new = line.split()
            merged[int(x)] = float(new)
        elif line.startswith("FALLBACK"):
            fallback = 1
        elif line.startswith("STATS"):
            stats = dict(kv.split("=") for kv in line.split()[1:])

    full = _W["full"]
    mism = 0
    for x in merged.keys() | full.keys():
        if merged.get(x, 0.0) != full.get(x, 0.0):
            mism += 1
    return dict(u=u, v=v, mismatches=mism, fallback=fallback,
                region=int(stats.get("region", -1)),
                pinned=int(stats.get("pinned", -1)),
                rounds=int(stats.get("rounds", -1)),
                pops=int(stats.get("pops", -1)),
                changed=int(stats.get("changed", -1)),
                insert_us=float(stats.get("insert_us", -1)))


def v3_static_time_us(graph_path, s, workdir):
    """Full-recompute baseline: V3 build+peel core time in us."""
    env = os.environ.copy()
    env.update({"PIVOTER_RUN_ST_V3": "1", "OMP_NUM_THREADS": "1"})
    out = subprocess.run(["./build/bin/degeneracy_cliques", graph_path, "1", str(s)],
                         env=env, capture_output=True, text=True, timeout=1800)
    txt = out.stdout + out.stderr
    build = re.search(r"ST_V3 Build took:\s*([\d.]+)\s*ms", txt)
    peel = re.search(r"STV3_PEEL_US:\s*(\d+)", txt)
    return float(build.group(1)) * 1000.0 + float(peel.group(1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--edges", type=int, default=300)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="dynamic_locality_out")
    ap.add_argument("--tmpdir", default=None)
    args = ap.parse_args()

    gname = Path(args.graph).stem
    workdir = args.tmpdir or tempfile.mkdtemp(prefix="dynins_")
    Path(workdir).mkdir(parents=True, exist_ok=True)
    Path(args.out).mkdir(parents=True, exist_ok=True)

    n, data_lines = read_edge_file(args.graph)
    print(f"[{gname} s={args.s}] n={n}, m={len(data_lines)}", flush=True)
    full_cores = run_v3(args.graph, args.s, "full", workdir)
    static_us = v3_static_time_us(args.graph, args.s, workdir)
    print(f"  static V3 build+peel = {static_us/1000.0:.1f} ms", flush=True)

    rng = random.Random(args.seed)
    picks = rng.sample(range(len(data_lines)), args.edges)

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers,
                             initializer=_init_worker,
                             initargs=(args.graph, args.s, workdir, full_cores)) as ex:
        for i, r in enumerate(ex.map(one_edge, picks)):
            rows.append(r)
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(picks)}", flush=True)

    csv_path = Path(args.out) / f"insert_{gname}_s{args.s}.csv"
    keys = ["u", "v", "mismatches", "fallback", "region", "pinned", "rounds",
            "pops", "changed", "insert_us"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys + ["error"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {csv_path}")

    ok = [r for r in rows if "error" not in r]
    errs = len(rows) - len(ok)
    fallbacks = sum(1 for r in ok if r.get("fallback"))
    ok = [r for r in ok if not r.get("fallback")]
    bad = sum(1 for r in ok if r["mismatches"] > 0)
    times = sorted(r["insert_us"] for r in ok)
    regions = sorted(r["region"] for r in ok)
    rounds = sorted(r["rounds"] for r in ok)
    def pct(a, q): return a[min(len(a) - 1, int(q * len(a)))]
    print(f"\n===== INSERT PROTOTYPE  {gname} s={args.s}  ({len(ok)} edges) =====")
    print(f"  errors: {errs}   fallbacks: {fallbacks}   "
          f"VERIFY FAILURES: {bad}  (must be 0)")
    print(f"  insert time us: median={pct(times,0.5):.0f}  p90={pct(times,0.9):.0f} "
          f" p99={pct(times,0.99):.0f}  max={times[-1]:.0f}")
    print(f"  region size: median={pct(regions,0.5)}  p90={pct(regions,0.9)} "
          f" max={regions[-1]}")
    print(f"  expand rounds: median={pct(rounds,0.5)}  max={rounds[-1]}")
    print(f"  static V3 recompute = {static_us:.0f} us")
    print(f"  SPEEDUP vs full recompute: median={static_us/max(1,pct(times,0.5)):.0f}x "
          f" p90={static_us/max(1,pct(times,0.9)):.0f}x  worst={static_us/max(1,times[-1]):.0f}x")


if __name__ == "__main__":
    main()
