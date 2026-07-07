#!/usr/bin/env python3
# §140 E4/E5: sample k random r-cliques from a graph (greedy common-neighborhood extension)
# for the nsi_query latency bench. Uniformity is irrelevant for latency; validity matters.
#   python3 gen_queries.py <graph.edges> <r> <k> <out.txt>
import sys, random

graph, r, k, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
adj = {}
with open(graph) as f:
    f.readline()
    for ln in f:
        p = ln.split()
        if len(p) < 2: continue
        a, b = int(p[0]), int(p[1])
        if a == b: continue
        adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
vs = [v for v in adj if len(adj[v]) >= r - 1]
random.seed(20260708)
got = 0
with open(out, "w") as f:
    tries = 0
    while got < k and tries < k * 200:
        tries += 1
        v = random.choice(vs)
        clique = [v]; cand = set(adj[v])
        while len(clique) < r and cand:
            u = random.choice(tuple(cand))
            clique.append(u); cand &= adj[u]
        if len(clique) == r:
            f.write(" ".join(map(str, clique)) + "\n"); got += 1
print(f"sampled {got}/{k} r-cliques -> {out}")
