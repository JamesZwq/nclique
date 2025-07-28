#!/usr/bin/env python3
"""
Build Succinct Clique Tree (SCT) exactly per Jain et al.
LaTeX pseudocode  ─  no efficiency tricks, pure Python.

Usage:
    python build_sct.py graph.txt
"""

import sys, heapq
# ------------------------------------------------------------
# 1. read edge list
# ------------------------------------------------------------
def read_graph(path):
    with open(path, "r") as f:
        n, m = map(int, f.readline().split())
        G = {i: set() for i in range(n)}
        for line in f:
            if not line.strip():
                continue
            u, v = map(int, line.split())
            G[u].add(v)
            G[v].add(u)
    return G

# ------------------------------------------------------------
# 2. degeneracy ordering (naïve heap O(m log n))
# ------------------------------------------------------------
def degeneracy_order(G):
    remaining = set(G)
    deg = {v: len(G[v]) for v in G}
    heap = [(deg[v], v) for v in G]
    heapq.heapify(heap)
    order = []
    while remaining:
        while True:
            d, v = heapq.heappop(heap)
            if v in remaining and d == deg[v]:
                break
        order.append(v)
        remaining.remove(v)
        for nbr in G[v]:
            if nbr in remaining:
                deg[nbr] -= 1
                heapq.heappush(heap, (deg[nbr], nbr))
    return order                      # v₁ … vₙ

# ------------------------------------------------------------
# 3. build SCT paths
# ------------------------------------------------------------
def build_sct(G):
    order = degeneracy_order(G)
    pos   = {v: i for i, v in enumerate(order)}
    # orient edges: N⁺(v) = {w : pos[w] > pos[v]}
    Nplus = {v: {w for w in G[v] if pos[w] > pos[v]} for v in G}

    paths = []

    def expand(H, P, Pstar):
        if not P:                          # leaf: output path ⟨H, P⋆⟩
            paths.append((tuple(sorted(H)), tuple(sorted(Pstar))))
            return
        # pivot u ∈ P maximising |P ∩ N(u)|
        u = max(P, key=lambda x: len(P & Nplus[x]))
        # iterate non-neighbours of pivot
        for v in list(P - Nplus[u]):
            P.remove(v)                    # delete v at this depth
            if v == u:                     # pivot branch → label v as pivot
                expand(H, P & Nplus[v], Pstar | {v})
            else:                          # hold branch  → label v as hold
                expand(H | {v}, P & Nplus[v], Pstar)

    for v in order:                        # one root per vertex
        expand({v}, set(Nplus[v]), set())
    return paths

# ------------------------------------------------------------
# 4. driver
# ------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python build_sct.py graph.txt")
        sys.exit(1)

    G = read_graph(sys.argv[1])
    paths = build_sct(G)

    print(f"Total paths: {len(paths)}\n")
    for H, P in paths:
        print(f"H={H}  P={P}")