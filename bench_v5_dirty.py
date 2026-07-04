#!/usr/bin/env python3
"""
v5 adjacency-guided dirty-set pricing — the fair last test for v5.

Window experiment killed the O-position walk (R* scattered in the order).
The surviving question: is the CLIQUE-ADJACENCY-guided repair local?  Its
cost is dominated by the vertices it must examine = R* plus the clique-boundary
(neighbors it checks and rejects).  Compare that boundary to v4's discovery
flood (tested 30k-41k / edge).

Per change edge (from the locality detail CSVs; R* = changed vertices):
  raw_bnd   = |N_G'(R*) \\ R*|                (graph-neighbour boundary)
  clq_bnd   = |{ y in N_G'(R*)\\R* : exists x in R* with |N(x) cap N(y)| >= s-2 }|
              (co-clique-able boundary = vertices the adjacency walk actually
               examines; a plausible dirty-frontier)
  reported vs |R*| and vs v4's tested-flood.

If clq_bnd is ~ |R*| (small), the adjacency walk examines little => v5's
"order/adjacency-guided discovery" beats v4's flood, v5 has a real edge.
If clq_bnd blows up (hub neighbours), v5's walk is no better than v4.

Usage:
  python3 bench_v5_dirty.py --graph graphs/com-dblp.edges \\
      --detail dynamic_locality_out/detail_com-dblp_s5.csv --s 5
"""
from __future__ import annotations
import argparse, csv
from collections import defaultdict


def load_adj(path):
    adj = defaultdict(set)
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p) < 2:
                continue
            a, b = int(p[0]), int(p[1])
            if a != b:
                adj[a].add(b); adj[b].add(a)
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--detail", required=True)
    ap.add_argument("--s", type=int, required=True)
    args = ap.parse_args()
    adj = load_adj(args.graph)
    k = args.s - 2   # co-clique threshold: co-inhabit an s-clique => >= s-2 common

    per_edge = defaultdict(list)
    with open(args.detail) as f:
        for row in csv.DictReader(f):
            per_edge[(int(row["u"]), int(row["v"]))].append(int(row["x"]))

    rows = []
    for (u, v), rstar in per_edge.items():
        R = set(rstar)
        # G' adjacency: add the inserted edge
        adju = adj[u] | {v}; adjv = adj[v] | {u}
        def nbrs(x):
            if x == u:
                return adju
            if x == v:
                return adjv
            return adj[x]
        # raw neighbour boundary
        bnd = set()
        for x in R:
            bnd |= nbrs(x)
        bnd -= R
        # co-clique boundary: y shares >= k common neighbours with some x in R*
        clq = 0
        for y in bnd:
            ny = nbrs(y)
            hit = False
            for x in R:
                if x not in ny:
                    continue
                nx = nbrs(x)
                # |N(x) cap N(y)| >= k ?
                a, b = (nx, ny) if len(nx) <= len(ny) else (ny, nx)
                c = 0
                for w in a:
                    if w in b:
                        c += 1
                        if c >= k:
                            hit = True
                            break
                if hit:
                    break
            if hit:
                clq += 1
        rows.append((len(R), len(bnd), clq))

    def stat(i):
        a = sorted(r[i] for r in rows)
        def p(q): return a[min(len(a)-1, int(q*len(a)))]
        return f"med={p(.5)} p90={p(.9)} p99={p(.99)} max={a[-1]}"
    print(f"===== V5 DIRTY  {args.detail}  ({len(rows)} change edges, s={args.s}) =====")
    print(f"  |R*|:           {stat(0)}")
    print(f"  raw boundary:   {stat(1)}")
    print(f"  CLIQUE boundary (walk-examined frontier): {stat(2)}")
    ratios = sorted((r[2] / max(1, r[0])) for r in rows)
    def pr(q): return ratios[min(len(ratios)-1, int(q*len(ratios)))]
    print(f"  clq_bnd / |R*|: med={pr(.5):.1f} p90={pr(.9):.1f} max={ratios[-1]:.1f}")
    walk = sorted(r[0] + r[2] for r in rows)
    def pw(q): return walk[min(len(walk)-1, int(q*len(walk)))]
    print(f"  walk size (|R*|+clq_bnd): med={pw(.5)} p90={pw(.9)} p99={pw(.99)} max={walk[-1]}")


if __name__ == "__main__":
    main()
