#!/usr/bin/env python3
"""
Correctness + performance driver for the dynamic (1,s)-core DELETE prototype.

Per sampled edge e=(u,v) of G:
  cores(G) computed once (full).  Prototype gets G (contains e) + cores(G),
  deletes e, emits CHANGED (x, old, new) with new < old.
  merged = cores(G) patched with CHANGED.
  VERIFY: merged == V3(G-e) per-vertex (missing id = core 0).

Baseline: V3 static build+peel and peel-only on G-e (recompute-from-scratch).

Usage:
  python3 bench_dynamic_delete.py --graph graphs/com-dblp.edges --s 5 \
      --edges 300 --workers 6 --tmpdir <scratch>
"""
from __future__ import annotations
import argparse, csv, os, random, re, subprocess, tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from bench_dynamic_locality import read_edge_file, run_v3

PROTO = os.environ.get("DYN1S_PROTO", "./build/bin/dynamic_1s_core_del")
_W = {}


def _init(graph_path, s, workdir, full_cores):
    n, lines = read_edge_file(graph_path)
    _W.update(n=n, lines=lines, s=s, workdir=workdir, full=full_cores,
              gpath=graph_path)


def one_edge(idx):
    lines, n, s, wd = _W["lines"], _W["n"], _W["s"], _W["workdir"]
    p = lines[idx].split(); u, v = int(p[0]), int(p[1])
    tag = f"e{idx}_{os.getpid()}"
    # reference: V3(G - e)
    gme = os.path.join(wd, f"gme_{tag}.edges")
    with open(gme, "w") as f:
        f.write(f"{n} {len(lines)-1}\n")
        for i, ln in enumerate(lines):
            if i != idx:
                f.write(ln + "\n")
    try:
        ref = run_v3(gme, s, tag + "r", wd)
    finally:
        os.unlink(gme)
    # prototype input: full G (contains e) + cores(G)
    cf = os.path.join(wd, f"cg_{tag}.tsv")
    with open(cf, "w") as f:
        for x, c in _W["full"].items():
            f.write(f"{x}\t{c:.0f}\n")
    try:
        proc = subprocess.run([PROTO, _W["gpath"], str(s), cf, str(u), str(v)],
                              capture_output=True, text=True, timeout=600)
    finally:
        os.unlink(cf)
    if proc.returncode != 0:
        return dict(u=u, v=v, error=proc.stderr[-200:])
    merged = dict(_W["full"])
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
    mism = 0
    for x in merged.keys() | ref.keys():
        if merged.get(x, 0.0) != ref.get(x, 0.0):
            mism += 1
    return dict(u=u, v=v, mismatches=mism, fallback=fallback,
                region=int(stats.get("region", -1)),
                rounds=int(stats.get("rounds", -1)),
                changed=int(stats.get("changed", -1)),
                insert_us=float(stats.get("insert_us", -1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--edges", type=int, default=300)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tmpdir", default=None)
    args = ap.parse_args()
    gname = Path(args.graph).stem
    wd = args.tmpdir or tempfile.mkdtemp(prefix="del_")
    Path(wd).mkdir(parents=True, exist_ok=True)
    n, lines = read_edge_file(args.graph)
    print(f"[{gname} s={args.s}] n={n} m={len(lines)}", flush=True)
    full = run_v3(args.graph, args.s, "full", wd)
    rng = random.Random(args.seed)
    picks = rng.sample(range(len(lines)), args.edges)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init,
                             initargs=(args.graph, args.s, wd, full)) as ex:
        for i, r in enumerate(ex.map(one_edge, picks)):
            rows.append(r)
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(picks)}", flush=True)
    ok = [r for r in rows if "error" not in r]
    errs = len(rows) - len(ok)
    fb = sum(1 for r in ok if r.get("fallback"))
    ok2 = [r for r in ok if not r.get("fallback")]
    bad = sum(1 for r in ok2 if r["mismatches"] > 0)
    t = sorted(r["insert_us"] for r in ok2)
    reg = sorted(r["region"] for r in ok2)
    def pct(a, q): return a[min(len(a)-1, int(q*len(a)))] if a else 0
    print(f"\n===== DELETE PROTOTYPE  {gname} s={args.s}  ({len(ok2)} edges) =====")
    print(f"  errors: {errs}   fallbacks: {fb}   VERIFY FAILURES: {bad}  (must be 0)")
    print(f"  delete time us: median={pct(t,.5):.0f} p90={pct(t,.9):.0f} "
          f"p99={pct(t,.99):.0f} max={t[-1] if t else 0:.0f}")
    print(f"  region: median={pct(reg,.5)} p90={pct(reg,.9)} max={reg[-1] if reg else 0}")


if __name__ == "__main__":
    main()
