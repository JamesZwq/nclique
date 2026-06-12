#!/usr/bin/env python3
"""
verify_hierarchy_brute.py — brute-force oracle for the vertex-level
(r,s)-nucleus hierarchy emitted by PIVOTER_RUN_REGION_V3LM_HIER.

Ground truth per level k:
  1. kappa per r-clique via the standard peel (witness s-cliques die
     when any constituent r-clique is peeled).
  2. Level-k witnesses: s-cliques S with min_{R subset S} kappa(R) >= k.
  3. Nuclei: connected components of witnesses under >=r vertex sharing.
  4. Vertex-level components: union nuclei whose vertex sets intersect
     (matches the class-DSU output: a shared vertex means a shared class).
Compare sorted component sizes and total vertices per printed level.

Usage: python3 scripts/verify_hierarchy_brute.py [trials]
"""
import itertools, random, re, subprocess, sys, tempfile, os

BIN = "./build/bin/degeneracy_cliques"


def cliques_of_size(adj, n, k):
    out = []
    def extend(clq, cand):
        if len(clq) == k:
            out.append(tuple(clq)); return
        for v in list(cand):
            if all(v in adj[u] for u in clq):
                extend(clq + [v], [w for w in cand if w > v])
    extend([], list(range(n)))
    return out


def brute_kappa(adj, n, r, s):
    rcl = cliques_of_size(adj, n, r)
    scl = cliques_of_size(adj, n, s)
    rid = {c: i for i, c in enumerate(rcl)}
    subs = [[rid[t] for t in itertools.combinations(S, r)] for S in scl]
    sup = [0] * len(rcl)
    for ss in subs:
        for t in ss: sup[t] += 1
    alive_r = [True] * len(rcl)
    alive_s = [True] * len(scl)
    kappa = [0] * len(rcl)
    k = 0
    remaining = len(rcl)
    while remaining:
        t = min((i for i in range(len(rcl)) if alive_r[i]), key=lambda i: sup[i])
        k = max(k, sup[t])
        kappa[t] = k
        alive_r[t] = False
        remaining -= 1
        for si, ss in enumerate(subs):
            if alive_s[si] and t in ss:
                alive_s[si] = False
                for o in ss:
                    if alive_r[o]: sup[o] -= 1
    return rcl, scl, kappa, rid


def vertex_components(adj, n, r, s, rcl, scl, kappa, rid, level):
    wit = [S for S in scl
           if min(kappa[rid[t]] for t in itertools.combinations(S, r)) >= level]
    if not wit:
        return []
    parent = list(range(len(wit)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i in range(len(wit)):
        for j in range(i + 1, len(wit)):
            if len(set(wit[i]) & set(wit[j])) >= r:
                parent[find(i)] = find(j)
    nuclei = {}
    for i, S in enumerate(wit):
        nuclei.setdefault(find(i), set()).update(S)
    vsets = list(nuclei.values())
    # merge vertex-overlapping nuclei
    merged = True
    while merged:
        merged = False
        for i in range(len(vsets)):
            for j in range(i + 1, len(vsets)):
                if vsets[i] & vsets[j]:
                    vsets[i] |= vsets[j]; del vsets[j]; merged = True; break
            if merged: break
    return sorted((len(v) for v in vsets), reverse=True)


def run_binary(path, r, s):
    env = dict(os.environ, PIVOTER_RUN_REGION_V3LM_HIER="1")
    p = subprocess.run([BIN, path, str(r), str(s), "degen"],
                       env=env, capture_output=True, text=True, timeout=120)
    out = {}
    for m in re.finditer(r"core>=(\d+)\s+verts=(\d+)\s+comps=(\d+)\s+sizes=\[([\d,\.]*)\]",
                         p.stdout):
        lvl = int(m.group(1))
        sizes = [int(x) for x in m.group(4).split(",") if x and x != "..."]
        out[lvl] = sizes
    return out


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    random.seed(20260613)
    fails = 0
    for t in range(trials):
        n = random.randint(9, 12)
        p = random.choice([0.5, 0.65, 0.8])
        r, s = random.choice([(3, 4), (3, 5)])
        adj = {v: set() for v in range(n)}
        edges = []
        for a, b in itertools.combinations(range(n), 2):
            if random.random() < p:
                adj[a].add(b); adj[b].add(a); edges.append((a, b))
        with tempfile.NamedTemporaryFile("w", suffix=".edges", delete=False) as f:
            f.write(f"{n} {len(edges)}\n")
            for a, b in edges: f.write(f"{a} {b}\n")
            path = f.name
        rcl, scl, kappa, rid = brute_kappa(adj, n, r, s)
        got = run_binary(path, r, s)
        ok = True
        for lvl, sizes in got.items():
            want = vertex_components(adj, n, r, s, rcl, scl, kappa, rid, lvl)
            if sizes != want:
                ok = False
                print(f"FAIL t={t} n={n} p={p} (r,s)=({r},{s}) level={lvl}: "
                      f"binary={sizes} brute={want}  [{path}]")
        if ok:
            os.unlink(path)
            print(f"ok t={t} n={n} p={p} ({r},{s}) levels={sorted(got)}")
        else:
            fails += 1
    print(f"\n{trials - fails}/{trials} passed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
