#!/usr/bin/env python3
"""
probe_false_twins.py — census of true-twin vs false-twin classes.

True twins : N[u] = N[v]  (closed; adjacent pairs)  — what profile
             classes can capture today (profile classes are coarser
             in principle but true twins are the vertex-level core).
False twins: N(u) = N(v), u !~ v  (open; non-adjacent pairs) — kappa-
             preserving via automorphism but invisible to profile
             classes (the cocktail-party gap).

Reports, per graph: #vertices, #nontrivial classes of each kind,
#vertices removable by quotienting (sum |class|-1), both globally and
restricted to vertices with degree >= s-1 (the only vertices that can
appear in any s-clique; default s=4 proxy).

Usage: python3 scripts/probe_false_twins.py [s] graph1.edges graph2.edges ...
"""
import sys
from collections import defaultdict


def census(path, s):
    adj = defaultdict(set)
    with open(path) as f:
        first = True
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                u, v = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            if first:
                first = False
                # many graphs carry a "n m" header line; detect: if this
                # pair never reappears as an edge it is harmless anyway
            if u == v:
                continue
            adj[u].add(v)
            adj[v].add(u)

    n = len(adj)
    deg_min = s - 1

    def group(closed):
        groups = defaultdict(list)
        for v, nb in adj.items():
            key = frozenset(nb | {v}) if closed else frozenset(nb)
            groups[hash(key)].append(v)
        return groups

    out = {}
    for kind, closed in (("true", True), ("false", False)):
        groups = group(closed)
        nontriv = removable = 0
        nontriv_hi = removable_hi = 0
        largest = 0
        for _, vs in groups.items():
            if len(vs) < 2:
                continue
            # verify (hash collision guard + for open: exclude pairs that
            # are adjacent, i.e. true-twin contamination of the open key)
            ref = adj[vs[0]] if not closed else (adj[vs[0]] | {vs[0]})
            ok = [v for v in vs
                  if (adj[v] if not closed else (adj[v] | {v})) == ref]
            if len(ok) < 2:
                continue
            nontriv += 1
            removable += len(ok) - 1
            largest = max(largest, len(ok))
            hi = [v for v in ok if len(adj[v]) >= deg_min]
            if len(hi) >= 2:
                nontriv_hi += 1
                removable_hi += len(hi) - 1
        out[kind] = (nontriv, removable, nontriv_hi, removable_hi, largest)
    return n, out


def main():
    args = sys.argv[1:]
    s = 4
    if args and args[0].isdigit():
        s = int(args[0])
        args = args[1:]
    print(f"(degree filter: deg >= s-1 = {s-1})")
    hdr = (f"{'graph':<16}{'n':>8} | "
           f"{'trueCls':>8}{'rmv':>9}{'rmv(hi)':>9}{'max':>6} | "
           f"{'falseCls':>9}{'rmv':>9}{'rmv(hi)':>9}{'max':>6}")
    print(hdr)
    for path in args:
        name = path.split("/")[-1].replace(".edges", "")
        n, out = census(path, s)
        t = out["true"]; fa = out["false"]
        print(f"{name:<16}{n:>8} | "
              f"{t[0]:>8}{t[1]:>9}{t[3]:>9}{t[4]:>6} | "
              f"{fa[0]:>9}{fa[1]:>9}{fa[3]:>9}{fa[4]:>6}", flush=True)


if __name__ == "__main__":
    main()
