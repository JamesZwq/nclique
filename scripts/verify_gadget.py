import itertools
# G = K_{8,8} (left L0..L7, right R0..R7) DISJOINT K_6 (w0..w5), plus u=L0 joined to w0,w1,w2
d, m = 8, 6
L = [f"L{i}" for i in range(d)]; R = [f"R{i}" for i in range(d)]; W = [f"w{i}" for i in range(m)]
V = L + R + W
E = set()
for a in L:
    for b in R: E.add(frozenset((a, b)))          # K_{8,8}
for a, b in itertools.combinations(W, 2): E.add(frozenset((a, b)))   # K_6
u = "L0"
for j in range(3): E.add(frozenset((u, f"w{j}")))  # u joined to w0,w1,w2
adj = {v: set() for v in V}
for e in E:
    a, b = tuple(e); adj[a].add(b); adj[b].add(a)
print(f"|V|={len(V)} |E|={len(E)} deg(u)={len(adj[u])} deg(w0)={len(adj['w0'])} deg(w3)={len(adj['w3'])}")

def support(v, alive, s):
    """# of s-cliques containing v, among alive vertices"""
    nb = [x for x in adj[v] if x in alive]
    cnt = 0
    for comb in itertools.combinations(nb, s - 1):
        if all(frozenset((a, b)) in E for a, b in itertools.combinations(comb, 2)): cnt += 1
    return cnt

def peel(s):
    """exact (1,s)-nucleus core numbers by peeling"""
    alive = set(V); core = {}; level = 0
    while alive:
        sup = {v: support(v, alive, s) for v in alive}
        v = min(alive, key=lambda x: sup[x])
        level = max(level, sup[v]); core[v] = level; alive.discard(v)
    return core

c2 = peel(2); c3 = peel(3)
print(f"\nkappa_(1,2): u={c2[u]}  w0={c2['w0']} w1={c2['w1']} w2={c2['w2']}  w3={c2['w3']}  L1={c2['L1']} R0={c2['R0']}")
print(f"kappa_(1,3): u={c3[u]}  w0={c3['w0']} w1={c3['w1']} w2={c3['w2']}  w3={c3['w3']}  L1={c3['L1']} R0={c3['R0']}")
print(f"\nFABLE CLAIMED: k12(u)=8 k12(w)=5 | k13(u)=3 k13(w)=10")
print(f"REVERSAL in cell(1,2): u({c2[u]}) > w0({c2['w0']}) ? {c2[u] > c2['w0']}")
print(f"REVERSAL in cell(1,3): u({c3[u]}) < w0({c3['w0']}) ? {c3[u] < c3['w0']}")
print(f"=> INTERACTING ORDER REVERSAL EXISTS: {(c2[u] > c2['w0']) and (c3[u] < c3['w0'])}")
# interaction check
tri_u_w = [(u,'w%d'%i,'w%d'%j) for i,j in itertools.combinations(range(3),2)]
print(f"u-w triangles (interaction witnesses in cell 1,3): {len(tri_u_w)}  edges u-w (cell 1,2): 3")
