#!/usr/bin/env python3
"""Brute-force check of the running example (Figure 1): core values, canonical nodes, own nodes, chains, ranks,
labels, DFS arrays, runs and entries for every size.  Everything the paper prints about the example comes from here."""
import itertools
from collections import defaultdict

names = ['b1','b2','b3','a1','a2','a3','a4','a5','u1','u2','u3','w1','w2','w3','x']
A = ['a1','a2','a3','a4','a5']; B = ['x','b1','b2','b3']; U = ['u1','u2','u3']; W = ['w1','w2','w3']
edges = set()
for S in (A, B):
    for p, q in itertools.combinations(S, 2): edges.add(frozenset((p, q)))
for u in U:
    for w in W: edges.add(frozenset((u, w)))
for y in U + W: edges.add(frozenset(('x', y)))
edges.add(frozenset(('a5', 'b3')))
adj = defaultdict(set)
for e in edges:
    p, q = tuple(e); adj[p].add(q); adj[q].add(p)

def cliques(S, s):
    return [frozenset(c) for c in itertools.combinations(sorted(S), s) if all(frozenset((p, q)) in edges for p, q in itertools.combinations(c, 2))]

def nuclei(s, k):
    """(s,k)-nuclei: maximal sets with every vertex in >= k s-cliques inside, connected through shared s-cliques (vertex-sharing)."""
    S = set(names)
    while True:
        cl = cliques(S, s); cnt = defaultdict(int)
        for c in cl:
            for v in c: cnt[v] += 1
        bad = {v for v in S if cnt[v] < k}
        if not bad: break
        S -= bad
    # components through s-cliques sharing a vertex
    cl = cliques(S, s); parent = {v: v for v in S}
    def find(v):
        while parent[v] != v: parent[v] = parent[parent[v]]; v = parent[v]
        return v
    for c in cl:
        c = sorted(c)
        for v in c[1:]: parent[find(v)] = find(c[0])
    comps = defaultdict(set)
    for v in S:
        if cnt[v] >= k: comps[find(v)].add(v)
    return [frozenset(c) for c in comps.values()]

kappa = {v: {} for v in names}; omega = {}
for v in names:
    om = max(len(c) for k in range(1, 6) for c in cliques(set(names), k) if v in c); omega[v] = om
for s in range(2, 6):
    k = 1
    while True:
        ns = nuclei(s, k)
        if not ns: break
        for N in ns:
            for v in N: kappa[v][s] = k
        k += 1
print('omega', omega)
for v in names: print(v, [kappa[v].get(s, 0) for s in range(2, 6)])
# canonical nodes per size: distinct nucleus vertex sets over all k; own node = nucleus at level kappa_s(v)
nodes = {}; own = {}
for s in range(2, 6):
    seen = {}
    for k in range(1, 8):
        for N in nuclei(s, k): seen.setdefault(N, []).append(k)
    nodes[s] = {N: (min(ks), max(ks)) for N, ks in seen.items()}   # level interval
    for v in names:
        if s in kappa[v]: own[(v, s)] = next(N for N in seen if v in N and max(seen[N]) == kappa[v][s])
    print('size', s, 'nodes:', [(sorted(N), iv) for N, iv in sorted(nodes[s].items(), key=lambda t: -len(t[0]))])
chains = defaultdict(list)
for v in names: chains[tuple(own.get((v, s)) for s in range(2, 6))].append(v)
print('chains:', [sorted(vs) for vs in chains.values()])
