"""Independent brute-force check of the lemmas in
research/r1_skyline_index_20260918/THEORY.md on small random graphs.

Checks: F1, F2, F3, F6 (dominance => core and nucleus containment), F7
(skyline = maximal pairs), F8 (absorption), canonical intervals, L1, L2,
L3, and Theorem Q2 (the retrieval procedure returns exactly the nucleus,
each vertex once).  Exact integers throughout.  Python is fine here: the
graphs have at most 9 vertices.
"""
import itertools, random, sys
from math import comb

def C(a, b):
    return comb(a, b) if 0 <= b <= a else 0

def cascade(k, s):
    """s-cascade of k>=1: list of (a_j, j) with a_s > a_{s-1} > ... >= j"""
    out = []; rem = k; j = s
    while rem > 0:
        assert j >= 1
        a = j
        while C(a + 1, j) <= rem: a += 1
        out.append((a, j)); rem -= C(a, j); j -= 1
    return out

def shadow(k, s):
    if k == 0: return 0
    return sum(C(a, j - 1) for a, j in cascade(k, s))

# ---- shadow sanity (Kruskal-Katona is tight on colex initial segments)
def colex_check():
    for s in (2, 3, 4):
        universe = list(range(12))
        subsets = sorted(itertools.combinations(universe, s), key=lambda t: tuple(reversed(t)))
        for k in range(1, 61):
            fam = subsets[:k]
            sh = {frozenset(sub) for S in fam for sub in itertools.combinations(S, s - 1)}
            assert len(sh) == shadow(k, s), (s, k, len(sh), shadow(k, s))
    for s in range(2, 6):
        prev = 0
        for k in range(0, 200):
            v = shadow(k, s); assert v >= prev; prev = v
    assert shadow(3, 2) == 3 and shadow(4, 2) == 4 and shadow(21, 2) == 7 and shadow(35, 3) == 21 and shadow(21, 3) == 17 and shadow(22, 3) == 18
colex_check()

def random_graph(rng):
    n = rng.randint(4, 11)
    style = rng.random()
    edges = set()
    if style < 0.5:
        p = rng.uniform(0.35, 0.9)
        for u, v in itertools.combinations(range(n), 2):
            if rng.random() < p: edges.add((u, v))
    else:  # planted cliques joined sparsely
        parts = []; i = 0
        while i < n:
            m = n - i if n - i < 2 else rng.randint(2, min(5, n - i)); parts.append(list(range(i, i + m))); i += m
        for P in parts:
            for u, v in itertools.combinations(P, 2): edges.add((u, v))
        for _ in range(rng.randint(0, n)):
            u, v = rng.sample(range(n), 2); edges.add((min(u, v), max(u, v)))
    return n, edges

def analyse(n, edges):
    adj = {v: set() for v in range(n)}
    for u, v in edges: adj[u].add(v); adj[v].add(u)
    def is_clique(S): return all(v in adj[u] for u, v in itertools.combinations(S, 2))
    def cliques(s, W): return [frozenset(c) for c in itertools.combinations(sorted(W), s) if is_clique(c)]
    omega = {v: max([len(c) for c in (frozenset(x) for s in range(1, n + 1) for x in itertools.combinations(range(n), s)) if v in c and is_clique(c)]) for v in range(n)}
    S = max(omega.values())
    kappa = {}
    for s in range(2, S + 2):
        alive = set(range(n)); k = 0; val = {}
        while alive:
            cl = cliques(s, alive); supp = {v: 0 for v in alive}
            for c in cl:
                for v in c: supp[v] += 1
            v = min(alive, key=lambda x: (supp[x], x)); k = max(k, supp[v]); val[v] = k; alive.remove(v)
        kappa[s] = val
    def nuclei(s, k):
        W = {v for v in range(n) if kappa[s][v] >= k}
        par = {v: v for v in W}
        def f(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for c in cliques(s, W):
            c = sorted(c)
            for v in c[1:]: par[f(v)] = f(c[0])
        comps = {}
        for v in W: comps.setdefault(f(v), set()).add(v)
        return [frozenset(c) for c in comps.values()]
    return adj, omega, S, kappa, nuclei

def check(n, edges, stats):
    adj, omega, S, kappa, nuclei = analyse(n, edges)
    # F1, F2
    for v in range(n):
        for s in range(2, S + 1):
            if kappa[s + 1][v] >= 1:
                assert kappa[s][v] >= s, 'F1'
                assert kappa[s][v] >= shadow(kappa[s + 1][v], s), 'F2'
                stats['F2'] += 1
        # support interval
        pos = [s for s in range(2, S + 2) if kappa[s][v] > 0]
        assert pos == list(range(2, omega[v] + 1)) if omega[v] >= 2 else pos == []
    # canonical nodes per s
    T = {}
    for s in range(2, S + 1):
        levels = sorted({kappa[s][v] for v in range(n) if kappa[s][v] > 0}, reverse=True)
        nodes = {}  # vertexset -> set of k where it is a nucleus
        for k in range(1, (levels[0] if levels else 0) + 1):
            for N in nuclei(s, k):
                nodes.setdefault(N, set()).add(k)
        for N, ks in nodes.items():
            assert ks == set(range(min(ks), max(ks) + 1)), 'interval'
            assert max(ks) == min(kappa[s][v] for v in N), 'k_hi = min kappa'
        T[s] = {N: (min(ks), max(ks)) for N, ks in nodes.items()}
    def own(s, v):  # X_s(v)
        for N, (lo, hi) in T[s].items():
            if v in N and lo <= kappa[s][v] <= hi: return N
        raise AssertionError('no own node')
    def ahi(s, M):  # container in T_s of node M of T_{s+1}
        l = shadow(T[s + 1][M][1], s)
        cands = [N for N, (lo, hi) in T[s].items() if M <= N and lo <= l <= hi]
        assert len(cands) == 1, 'F3 container'
        return cands[0]
    # L1, L2, F3, F7, F8
    for v in range(n):
        if omega[v] < 2: continue
        sky = []
        for s in range(2, omega[v] + 1):
            assert T[s][own(s, v)][1] == kappa[s][v], 'L1'
            if s == omega[v]: sky.append(s); continue
            delta = kappa[s][v] - shadow(kappa[s + 1][v], s)
            if delta > 0: sky.append(s)
            else:
                assert own(s, v) == ahi(s, own(s + 1, v)), 'L2'; stats['L2'] += 1
        # F7: maximal pairs under dominance == skyline
        def dom(sp, kp, s, k):
            if sp < s: return False
            x = kp
            for t in range(sp - 1, s - 1, -1): x = shadow(x, t)
            return x >= k
        pairs = {s: kappa[s][v] for s in range(2, omega[v] + 1)}
        maximal = [s for s in pairs if not any(t != s and dom(t, pairs[t], s, pairs[s]) for t in pairs)]
        assert sorted(maximal) == sky, ('F7', maximal, sky)
        # F8
        eq = [s for s in range(2, omega[v] + 1) if kappa[s][v] == C(omega[v] - 1, s - 1)]
        if eq:
            assert eq == list(range(min(eq), omega[v] + 1)), 'F8'
            for s in range(min(eq), omega[v]):
                assert kappa[s][v] - shadow(kappa[s + 1][v], s) == 0, 'F8 delta'
        stats['vertices'] += 1
    # F6: dominance => containment of cores and nuclei (sample all occurring pairs)
    for s in range(2, S + 1):
        for sp in range(s, S + 1):
            for kp in sorted({kappa[sp][v] for v in range(n) if kappa[sp][v] > 0}):
                x = kp
                for t in range(sp - 1, s - 1, -1): x = shadow(x, t)
                for k in range(1, x + 1):
                    Wp = {v for v in range(n) if kappa[sp][v] >= kp}; W = {v for v in range(n) if kappa[s][v] >= k}
                    assert Wp <= W, 'F6 core'
                    for Np in nuclei(sp, kp):
                        assert sum(1 for N in nuclei(s, k) if Np <= N) == 1, 'F6 nucleus'
                    stats['F6'] += 1
    # Theorem Q2: simulate the index and the query
    sky_of = {}
    for v in range(n):
        if omega[v] < 2: continue
        sk = []
        for s in range(2, omega[v] + 1):
            if s == omega[v] or kappa[s][v] - shadow(kappa[s + 1][v], s) > 0: sk.append(s)
        sky_of[v] = sk
    buckets = {}  # (s, node) -> list of (v, gamma)
    for v, sk in sky_of.items():
        prev = 0
        for s in sk:
            buckets.setdefault((s, own(s, v)), []).append((v, prev)); prev = s
    def query(v, s, k):
        sk = sky_of[v]; sp = min(t for t in sk if t >= s)
        X = own(sp, v)
        for t in range(sp - 1, s - 1, -1): X = ahi(t, X)
        assert X == own(s, v), 'chain'
        N = X
        while True:  # walk up while parent's top >= k
            lo, hi = T[s][N]
            if lo <= k: break
            parent = [P for P, (plo, phi) in T[s].items() if N < P and phi == lo - 1]
            assert len(parent) == 1, 'parent'
            N = parent[0]
        adm = [M for M in T[s] if M <= N]  # N and descendants
        out = []
        t = s
        while adm and t <= S:
            for M in adm:
                for (u, g) in buckets.get((t, M), []):
                    if g < s: out.append(u)
            nxt = [M2 for M2 in T.get(t + 1, {}) if ahi(t, M2) in adm] if t + 1 in T else []
            adm = nxt; t += 1
        return N, out
    for v in range(n):
        if omega[v] < 2: continue
        for s in range(2, omega[v] + 1):
            for k in sorted({kappa[s][u] for u in range(n) if 0 < kappa[s][u] <= kappa[s][v]} | {1}):
                N, out = query(v, s, k)
                truth = next(M for M in nuclei(s, k) if v in M)
                assert N == truth, 'Q2 node'
                assert len(out) == len(set(out)) and set(out) == truth, ('Q2 output', sorted(out), sorted(truth))
                stats['Q2'] += 1
                # L3
                for u in range(n):
                    inside = kappa[s][u] >= 1 and own(s, u) <= N
                    assert inside == (u in truth), 'L3'

def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 20260918)
    stats = {'graphs': 0, 'vertices': 0, 'F2': 0, 'F6': 0, 'L2': 0, 'Q2': 0}
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    for _ in range(trials):
        n, edges = random_graph(rng)
        check(n, edges, stats); stats['graphs'] += 1
    # fixed examples: K_m plus pendant structures, the 14-vertex example
    ex = set()
    def cl(vs):
        for u, v in itertools.combinations(vs, 2): ex.add((min(u, v), max(u, v)))
    cl([0, 1, 2, 3, 4]); cl([4, 5, 6, 7]); cl([7, 8, 9])
    # octahedron omitted (n<=9 keeps cost low); add a bridge triangle
    cl([2, 8, 9])
    check(10, ex, stats); stats['graphs'] += 1
    print(stats)

main()
