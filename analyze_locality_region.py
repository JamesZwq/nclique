#!/usr/bin/env python3
"""
Affected-region analysis for dynamic (1,s)-core maintenance.

Input: detail CSV from bench_dynamic_locality.py (--detail rows:
u,v,x,old_core,new_core — one row per changed vertex per inserted edge).

For each inserted edge e=(u,v) we test the three ingredients of the
candidate soundness lemma:

  (1) HOP DISTANCE: BFS distance in G from {u,v} to each changed vertex.
      Small distances => the region is a small ball around the insertion.

  (2) INTERNAL CONNECTIVITY: in G[changed ∪ {u,v}], is every changed
      vertex connected to {u,v}?  The lemma predicts YES: every new
      survivor at level ℓ must share an s-clique with the new edge or
      with another new survivor, so the changed set ∪ {u,v} must be
      connected (s-clique connectivity implies plain connectivity in the
      induced subgraph).

  (3) LEVEL STRUCTURE: relation of old_core(x) to K = min(core_G-e'(u),
      core(v)) — for k-core the changed set is exactly {x : core = K};
      we measure how the analog behaves for s>=3 (expect: old cores lie
      in a band, not a single level).

Usage:
  python3 analyze_locality_region.py --graph graphs/com-dblp.edges \
      --detail dynamic_locality_out/detail_com-dblp_s5.csv
"""
from __future__ import annotations
import argparse, csv, sys
from collections import defaultdict, deque, Counter


def load_adj(path: str):
    adj = defaultdict(list)
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p) < 2:
                continue
            a, b = int(p[0]), int(p[1])
            if a == b:
                continue
            adj[a].append(b)
            adj[b].append(a)
    return adj


def bfs_dists(adj, sources, targets):
    """Hop distance from sources to each target (truncated when all found)."""
    want = set(targets)
    dist = {s: 0 for s in sources}
    found = {t: 0 for t in targets if t in dist}
    q = deque(sources)
    while q and len(found) < len(want):
        x = q.popleft()
        for y in adj[x]:
            if y not in dist:
                dist[y] = dist[x] + 1
                if y in want:
                    found[y] = dist[y]
                q.append(y)
    return [found.get(t, -1) for t in targets]  # -1 = unreachable


def connected_to_uv(adj, changed, u, v):
    """Within G[changed ∪ {u,v}], count changed vertices connected to {u,v}."""
    allowed = set(changed) | {u, v}
    seen = {u, v}
    q = deque([u, v])
    while q:
        x = q.popleft()
        for y in adj[x]:
            if y in allowed and y not in seen:
                seen.add(y)
                q.append(y)
    return sum(1 for x in changed if x in seen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--detail", required=True)
    args = ap.parse_args()

    adj = load_adj(args.graph)

    per_edge = defaultdict(list)  # (u,v) -> [(x, old, new)]
    with open(args.detail) as f:
        for row in csv.DictReader(f):
            per_edge[(int(row["u"]), int(row["v"]))].append(
                (int(row["x"]), int(row["old_core"]), int(row["new_core"])))

    dist_ctr = Counter()
    n_edges = len(per_edge)
    conn_ok_edges = 0
    conn_total_changed = 0
    conn_connected = 0
    band_stats = []

    for (u, v), changed in per_edge.items():
        xs = [x for (x, _, _) in changed if x != u and x != v]
        if xs:
            for d in bfs_dists(adj, [u, v], xs):
                dist_ctr[d] += 1
        allx = [x for (x, _, _) in changed]
        got = connected_to_uv(adj, allx, u, v)
        conn_total_changed += len(allx)
        conn_connected += got
        if got == len(allx):
            conn_ok_edges += 1
        # level band: old cores of changed relative to their own min/max
        olds = [o for (_, o, _) in changed]
        if olds:
            band_stats.append((min(olds), max(olds)))

    total = sum(dist_ctr.values())
    print(f"===== REGION ANALYSIS  {args.detail}  ({n_edges} edges w/ changes) =====")
    print(f"  hop distance from {{u,v}} of changed vertices (excl. u,v):")
    for d in sorted(dist_ctr):
        c = dist_ctr[d]
        label = "unreachable" if d < 0 else f"d={d}"
        print(f"    {label}: {c}  ({100.0*c/total:.1f}%)")
    print(f"  induced-connectivity check (lemma prediction: 100%):")
    print(f"    edges where ALL changed are connected to {{u,v}} in "
          f"G[changed∪{{u,v}}]: {conn_ok_edges}/{n_edges} "
          f"({100.0*conn_ok_edges/max(1,n_edges):.1f}%)")
    print(f"    changed vertices connected: {conn_connected}/{conn_total_changed} "
          f"({100.0*conn_connected/max(1,conn_total_changed):.1f}%)")
    if band_stats:
        widths = sorted(mx - mn for (mn, mx) in band_stats)
        print(f"  old-core band width per edge (max-min of changed old cores):")
        print(f"    median={widths[len(widths)//2]}  "
              f"p90={widths[int(0.9*len(widths))]}  max={widths[-1]}")


if __name__ == "__main__":
    main()
