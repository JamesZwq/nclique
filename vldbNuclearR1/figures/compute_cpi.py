"""
Replicate the CPI (Clique Path Index) construction as used by the SIGMOD
baseline (src/SDCT_Fused.hpp + src/misc.cpp::findBestPivot...).

Algorithm outline (matches C++):
  1. Process vertices in degeneracy order.
  2. For each vertex v, run recursive BK on its later-neighbors.
  3. At each recursion: Tomita pivot u = argmax_{w ∈ P ∪ X} |P ∩ N(w)|.
     Branching set = {u} ∪ (P - N(u)).
     For each candidate c in branching set:
       - if c == u: push to dropV (→ V_p)
       - else:      push to keepV (→ V_h)
  4. Emit leaf (keepV, dropV) when P becomes empty.

A leaf (V_h, V_p) represents all s-cliques V_h ∪ P' where P' ⊆ V_p of size
(s - |V_h|). Total count = C(|V_p|, s - |V_h|).
"""
from math import comb
from itertools import combinations
from collections import defaultdict

def build_graph():
    edges = [
        (1,2),(1,3),(1,4),(2,3),(2,4),(3,4),           # Block A (K4)
        (3,5),(3,6),(4,5),(4,6),(5,6),                  # Block B (K4)
        (6,7),(6,8),(7,8),                              # T1
        (2,9),(2,10),(9,10),                            # T2
    ]
    adj = {i: set() for i in range(1, 11)}
    for u, v in edges:
        adj[u].add(v); adj[v].add(u)
    return adj

def degeneracy_order(adj):
    cur_deg = {v: len(nbrs) for v, nbrs in adj.items()}
    alive = set(adj)
    order = []
    while alive:
        v = min(alive, key=lambda x: (cur_deg[x], x))
        order.append(v)
        alive.remove(v)
        for u in adj[v]:
            if u in alive:
                cur_deg[u] -= 1
    return order

def sdct_cpi(adj, s):
    order = degeneracy_order(adj)
    pos = {v: i for i, v in enumerate(order)}
    leaves = []

    def choose_pivot(P, X):
        # Matches C++ findBestPivotNonNeighborsDegeneracyCliques: pivot from P only.
        return max(P, key=lambda w: len(P & adj[w]))

    def bk(keepV, dropV, P, X):
        if not P:
            if keepV or dropV:
                leaves.append((tuple(keepV), tuple(dropV)))
            return
        u = choose_pivot(P, X)
        # branching set = (P - N(u)) ∪ {u}  (u always in P here; if u in X, add P - N(u) only)
        if u in P:
            branch = (P - adj[u]) | {u}
        else:
            branch = P - adj[u]
        # iterate deterministically for reproducibility
        P_work = set(P); X_work = set(X)
        for c in sorted(branch):
            if c not in P_work:
                continue
            newP = P_work & adj[c]
            newX = X_work & adj[c]
            if c == u:
                bk(keepV, dropV + [c], newP, newX)
            else:
                bk(keepV + [c], dropV, newP, newX)
            P_work.discard(c); X_work.add(c)

    # Outer degeneracy loop
    for i, v in enumerate(order):
        P_later = {u for u in adj[v] if pos[u] > i}
        X_earlier = {u for u in adj[v] if pos[u] < i}
        bk([v], [], P_later, X_earlier)

    # Filter: only keep leaves that can encode at least one s-clique.
    kept = [(h, p) for (h, p) in leaves
            if len(h) <= s and len(h) + len(p) >= s]
    return kept, order

def verify(adj, s, leaves):
    # (1) Total s-clique count should match brute-force.
    brute = sum(1 for combo in combinations(sorted(adj), s)
                if all(b in adj[a] for a in combo for b in combo if b != a))
    encoded = sum(comb(len(p), s - len(h)) for (h, p) in leaves)

    # (2) Per-vertex support:
    #   hold vertex v contributes C(|V_p|, s-|V_h|)
    #   pivot vertex v contributes C(|V_p|-1, (s-|V_h|)-1)
    supp = defaultdict(int)
    for h, p in leaves:
        eta = s - len(h)
        for v in h:
            supp[v] += comb(len(p), eta)
        for v in p:
            if eta >= 1:
                supp[v] += comb(len(p) - 1, eta - 1)

    # brute-force support
    brute_supp = defaultdict(int)
    for combo in combinations(sorted(adj), s):
        if all(b in adj[a] for a in combo for b in combo if b != a):
            for v in combo:
                brute_supp[v] += 1

    return brute, encoded, dict(supp), dict(brute_supp)

def simulate_peeling(adj, leaves, s, verbose=True, tie_break="pedagogical"):
    """
    tie_break:
      'smallest_id' — classical, pops smallest id first (→ Case A dominates)
      'pedagogical'  — prefers vertices that are pivot of some ALIVE leaf
                       (produces a trace that showcases Case A + C + C')
    """
    """
    Simulate the Immutable-CPI peeling algorithm and return the trace.
    Trace rows: (round, popped_vertex, core_value, events: list of path-update dicts)
    """
    # Initial per-path state
    rp = [len(p) for h, p in leaves]       # remaining pivots
    eta = [s - len(h) for h, p in leaves]  # needPivot
    alive = [rp[i] >= eta[i] for i in range(len(leaves))]
    # Per-vertex support from CPI
    supp = defaultdict(int)
    for i, (h, p) in enumerate(leaves):
        for v in h:
            supp[v] += comb(rp[i], eta[i])
        for v in p:
            if eta[i] >= 1:
                supp[v] += comb(rp[i] - 1, eta[i] - 1)
    # vertex -> list of (leaf_idx, role) where role ∈ {"hold", "pivot"}
    vtx_in = defaultdict(list)
    for i, (h, p) in enumerate(leaves):
        for v in h:
            vtx_in[v].append((i, "hold"))
        for v in p:
            vtx_in[v].append((i, "pivot"))

    peeled = set()
    core = {}
    trace = []
    round_no = 0
    cur_level = 0
    while len(peeled) < len(adj):
        # advance level to the min support among non-peeled with support ≥ cur_level
        remain = [v for v in adj if v not in peeled]
        if not remain: break
        min_supp = min(max(supp[v], 0) for v in remain)
        cur_level = max(cur_level, min_supp)
        # pop ONE vertex at current level (to see case-by-case trace)
        if tie_break == "pedagogical":
            def prefer(v):
                s_v = supp[v]
                # Prefer pivot of alive leaves first (to surface Case C)
                pivot_alive = any(alive[i] and role == "pivot" for (i, role) in vtx_in[v])
                return (s_v, 0 if pivot_alive else 1, v)
            victim = min(remain, key=prefer)
        else:
            victim = min(remain, key=lambda v: (supp[v], v))
        if supp[victim] > cur_level:
            # should not happen
            cur_level = supp[victim]
        round_no += 1
        core[victim] = cur_level
        peeled.add(victim)
        events = []
        for (i, role) in vtx_in[victim]:
            if not alive[i]:
                events.append({"leaf": i+1, "status": "skip (dead)"})
                continue
            old_rp = rp[i]
            new_rp = old_rp - (1 if role == "pivot" else 0)
            # Case A: hold removed OR new_rp < eta[i]  → path dies
            if role == "hold":
                # full delta to remaining (alive) pivots; all holds are "gone"
                # Δ_pivot^dagger = C(|V_p|-1, eta-1)  — but we actually peel by the
                # "remaining" pivot count
                Dp = comb(old_rp - 1, eta[i] - 1) if eta[i] >= 1 else 0
                # surviving holds (none in this leaf since hold is just one here,
                # but in general there could be multiple holds in same leaf):
                Dh = comb(old_rp, eta[i])
                # Apply deltas
                for v in leaves[i][0]:
                    if v != victim and v not in peeled:
                        supp[v] -= Dh
                for v in leaves[i][1]:
                    if v not in peeled:
                        supp[v] -= Dp
                alive[i] = False
                events.append({"leaf": i+1, "case": "A (hold removed)",
                               "Dh": Dh, "Dp": Dp, "path_dies": True})
            else:
                # pivot removed
                if new_rp < eta[i]:
                    # Case C': path dies
                    Dh_full = comb(old_rp, eta[i])
                    Dp_full = comb(old_rp - 1, eta[i] - 1) if eta[i] >= 1 else 0
                    for v in leaves[i][0]:
                        if v not in peeled:
                            supp[v] -= Dh_full
                    for v in leaves[i][1]:
                        if v != victim and v not in peeled:
                            supp[v] -= Dp_full
                    alive[i] = False
                    rp[i] = new_rp
                    events.append({"leaf": i+1, "case": "C' (pivot → death)",
                                   "p": f"{old_rp}→{new_rp}",
                                   "Dh": Dh_full, "Dp": Dp_full, "path_dies": True})
                else:
                    # Case C: pivot shrink, path alive
                    Dh = comb(old_rp, eta[i]) - comb(new_rp, eta[i])
                    Dp = (comb(old_rp - 1, eta[i] - 1) - comb(new_rp - 1, eta[i] - 1)) \
                         if eta[i] >= 1 else 0
                    for v in leaves[i][0]:
                        if v not in peeled:
                            supp[v] -= Dh
                    for v in leaves[i][1]:
                        if v != victim and v not in peeled:
                            supp[v] -= Dp
                    rp[i] = new_rp
                    events.append({"leaf": i+1, "case": "C (pivot shrink)",
                                   "p": f"{old_rp}→{new_rp}",
                                   "Dh": Dh, "Dp": Dp, "path_dies": False})
        trace.append({"round": round_no, "victim": victim, "core": cur_level,
                      "events": events, "supp_after": dict(supp)})
    return trace, core

if __name__ == "__main__":
    adj = build_graph()
    order = degeneracy_order(adj)
    print("Degeneracy order:", order)

    s = 3
    leaves, _ = sdct_cpi(adj, s)
    print(f"\nCPI for s={s}: {len(leaves)} leaves")
    for i, (h, p) in enumerate(leaves, 1):
        n = comb(len(p), s - len(h))
        print(f"  L{i}: V_h = {list(h)}, V_p = {list(p)}  →  C({len(p)},{s-len(h)}) = {n}")

    brute, encoded, supp, brute_supp = verify(adj, s, leaves)
    print(f"\nTotal encoded: {encoded}; brute-force count: {brute}; match: {encoded == brute}")

    print("\nInitial supports:")
    for v in sorted(adj):
        print(f"  v{v}: {supp.get(v, 0)}")

    trace, core = simulate_peeling(adj, leaves, s)
    print("\nPeeling trace:")
    for row in trace:
        print(f"\nRound {row['round']}: pop v{row['victim']}, core={row['core']}")
        for e in row["events"]:
            items = [f"{k}={v}" for k, v in e.items() if k != 'leaf']
            print(f"    L{e['leaf']}: " + ", ".join(items))

    print("\nFinal cores:")
    for v in sorted(adj):
        print(f"  v{v}: {core[v]}")
