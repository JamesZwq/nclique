#!/usr/bin/env python3
"""
edges_to_gbbs_adj.py — convert our .edges (SNAP-style edge list with an
"n m" header, 0-indexed, each undirected edge listed once) into the GBBS
AdjacencyGraph CSR text format that NucleusDecomposition_main reads.

GBBS AdjacencyGraph layout (symmetric/undirected):
    AdjacencyGraph
    <n>
    <m>                 # total directed endpoints = 2*|undirected edges kept|
    <offset[0]>         # n offset lines, offset[i] = start of vertex i in edges[]
    ...
    <offset[n-1]>
    <edge[0]>           # m target-vertex lines, neighbor lists sorted
    ...
    <edge[m-1]>

We must convert OUR exact .edges so the ARB comparison runs on the same
graph as the RegND* numbers; the stale .adj files in /data/wenqianz are a
different (same-named) dataset and must not be reused.

Usage: python3 edges_to_gbbs_adj.py in.edges out.adj
"""
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: edges_to_gbbs_adj.py <in.edges> <out.adj>")
    inp, outp = sys.argv[1], sys.argv[2]

    n = None
    # first pass: read header for n, collect edges
    adj = None
    with open(inp) as f:
        first = f.readline().split()
        # header is "n m"; if the first line is actually an edge (two ids
        # both < a plausible n) we still treat it as header only when it
        # has exactly 2 ints AND the file convention is header-present.
        n = int(first[0])
        adj = [set() for _ in range(n)]
        for line in f:
            p = line.split()
            if len(p) < 2:
                continue
            u, v = int(p[0]), int(p[1])
            if u == v:
                continue          # drop self-loops (GBBS dislikes them)
            if u >= n or v >= n:
                # vertex id beyond header n: grow (header was a lower bound)
                m = max(u, v)
                adj.extend(set() for _ in range(m + 1 - len(adj)))
                n = len(adj)
            adj[u].add(v)
            adj[v].add(u)         # symmetric

    offsets = [0] * n
    total = 0
    for i in range(n):
        offsets[i] = total
        total += len(adj[i])

    with open(outp, "w") as g:
        g.write("AdjacencyGraph\n")
        g.write(f"{n}\n")
        g.write(f"{total}\n")
        g.write("\n".join(map(str, offsets)))
        g.write("\n")
        # edges, neighbor lists sorted
        buf = []
        for i in range(n):
            if adj[i]:
                buf.append("\n".join(map(str, sorted(adj[i]))))
        g.write("\n".join(buf))
        g.write("\n")

    print(f"{inp} -> {outp}: n={n} m(directed)={total} "
          f"|E|={total // 2}", flush=True)


if __name__ == "__main__":
    main()
