#!/usr/bin/env python3
"""
v5 repair-walk locality — decisive kill/confirm.

The repair walk restarts the peel from the earliest bumped seed p and walks
the O-suffix, touching only "dirty" vertices.  Its cost is local IFF the
changed set R* sits in a short O-window near p (few non-changing vertices to
step through).  If R* is scattered across the whole suffix, an order-guided
walk degenerates and v5 cannot hit stable-1000x.

Per CHANGE edge e=(u,v):
  base = G-e; V3 dumps base core[] + pop-rank[] (a total order O); full = V3(G).
  R* = changed vertices.  Split:
    R*_zero    = base_core 0 (absent rank) risers  -- trivially local (isolated
                 vertices joining their first clique; rank ~ front)
    R*_counted = base_core > 0 changed vertices, each with a real O-rank.
  For R*_counted (when nonempty):
    span      = max_rank - min_rank        (spread of R* in the order)
    cnt       = |R*_counted|
    density   = cnt / (span+1)             (1.0 = perfectly contiguous)
    p_rank    = min_rank                    (walk start proxy)
    suffix    = N_counted - p_rank          (order length after p)
    tail_frac = span / suffix               (does R* reach the order's end?)

Verdict:
  span small / density high / tail_frac small  => walk is LOCAL, v5 viable
  span ~ suffix (tail_frac ~ 1)                => R* scattered, walk degenerates

Usage:
  python3 bench_v5_window.py --graph graphs/com-dblp.edges --s 5 \
      --edges 300 --workers 6 --tmpdir <scratch>
"""
from __future__ import annotations
import argparse, os, random, tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from bench_dynamic_locality import read_edge_file, run_v3
from bench_v5_certificate import run_v3_order

_W = {}


def _init(graph_path, s, workdir, full_cores):
    n, lines = read_edge_file(graph_path)
    _W.update(n=n, lines=lines, s=s, workdir=workdir, full=full_cores)


def one_edge(idx):
    lines, n, s, wd = _W["lines"], _W["n"], _W["s"], _W["workdir"]
    p = lines[idx].split(); u, v = int(p[0]), int(p[1])
    tag = f"e{idx}_{os.getpid()}"
    gf = os.path.join(wd, f"g_{tag}.edges")
    with open(gf, "w") as f:
        f.write(f"{n} {len(lines)-1}\n")
        for i, ln in enumerate(lines):
            if i != idx:
                f.write(ln + "\n")
    try:
        bcore, brank = run_v3_order(gf, s, tag, wd)
    finally:
        os.unlink(gf)
    full = _W["full"]
    rstar = [x for x in set(bcore) | set(full) if bcore.get(x, 0.0) != full.get(x, 0.0)]
    if not rstar:
        return dict(change=0)
    counted = [x for x in rstar if x in brank]           # base rank exists
    zero = [x for x in rstar if x not in brank]          # base core 0 risers
    N = len(brank)                                        # counted-order length
    res = dict(change=1, nrstar=len(rstar), nzero=len(zero), ncounted=len(counted))
    if counted:
        ranks = [brank[x] for x in counted]
        lo, hi = min(ranks), max(ranks)
        span = hi - lo
        suffix = N - lo
        res.update(span=span, cnt=len(counted),
                   density=len(counted) / (span + 1),
                   suffix=suffix, tail_frac=(span / suffix) if suffix else 0.0)
    return res


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
    wd = args.tmpdir or tempfile.mkdtemp(prefix="v5w_")
    Path(wd).mkdir(parents=True, exist_ok=True)
    n, lines = read_edge_file(args.graph)
    print(f"[{gname} s={args.s}] n={n}", flush=True)
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

    ch = [r for r in rows if r.get("change")]
    wc = [r for r in ch if "span" in r]            # have counted R*
    def stat(key, arr):
        a = sorted(r[key] for r in arr)
        if not a:
            return "n/a"
        def p(q): return a[min(len(a)-1, int(q*len(a)))]
        return f"med={p(.5):.3g} p90={p(.9):.3g} p99={p(.99):.3g} max={a[-1]:.3g}"
    print(f"\n===== V5 WINDOW  {gname} s={args.s}  ({len(ch)} change edges) =====")
    print(f"  change edges with counted R*: {len(wc)}  "
          f"(zero-only risers: {len(ch)-len(wc)})")
    print(f"  |R*|:        {stat('nrstar', ch)}")
    print(f"  |R*_counted|:{stat('cnt', wc)}")
    print(f"  O-SPAN of R*: {stat('span', wc)}")
    print(f"  density(cnt/span): {stat('density', wc)}")
    print(f"  suffix len (N-p): {stat('suffix', wc)}")
    print(f"  TAIL_FRAC (span/suffix): {stat('tail_frac', wc)}")
    # decisive rollup: fraction of change edges whose R* is 'tight'
    tight = sum(1 for r in wc if r["span"] <= 10 * r["cnt"])
    print(f"  R* tight (span <= 10*cnt): {tight}/{len(wc)} "
          f"({100.0*tight/max(1,len(wc)):.0f}%)")


if __name__ == "__main__":
    main()
