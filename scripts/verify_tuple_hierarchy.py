#!/usr/bin/env python3
"""Brute-force verifier for the TUPLE-NATIVE (r,s)-nucleus hierarchy emitted by
region_native_sct_peel.cpp (the a_Y engine) under PIVOTER_DUMP_HIER.

What it checks, on many tiny random graphs:

  1. CORE CROSS-CHECK: brute-force r-clique cores -> core-distribution (weighted
     by count) must equal a_Y's printed core distribution. (The hierarchy pass
     must not change the peel result.)

  2. HIERARCHY EQUIVALENCE: the tuple-native forest must represent the SAME
     nuclei as the ground-truth r-clique-level forest. We compare a canonical
     "nucleus profile": the multiset over every distinct core-level k of
     (k, component_size) for each connected component of the subgraph induced by
     {r-cliques with core >= k} where adjacency = share an s-clique. component
     size = #r-cliques in the component. This multiset encodes births, sizes and
     the full parent/child nesting of the merge forest. Two forests represent the
     same nuclei iff their profiles are identical multisets.

     - GT profile is computed directly from r-clique cores + adjacency.
     - a_Y profile is reconstructed from the dumped merge-tree CSV.
     We also self-check the reconstruction by building the GT merge tree the same
     way and confirming its reconstructed profile == the direct GT profile.

Run:  python3 scripts/verify_tuple_hierarchy.py [N_INSTANCES] [BINARY]
Exit code 0 and "VIOLATIONS: 0" on success.
"""
import itertools
import os
import random
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

BINARY = "region_native/region_native_sct_peel"


# ----------------------------- brute force -----------------------------------
def all_cliques_of_size(adj, n, k):
    """Yield every k-clique (sorted tuple) in the graph given by adj (set of sets)."""
    verts = range(n)
    for combo in itertools.combinations(verts, k):
        ok = True
        for i in range(k):
            ai = adj[combo[i]]
            for j in range(i + 1, k):
                if combo[j] not in ai:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            yield combo


def brute_force_nucleus(adj, n, r, s):
    """Return (rcliques, core) where core[rc] = (r,s)-nucleus core number.

    Support of an r-clique R = number of s-cliques containing R. Peel: repeatedly
    remove the r-clique of minimum current support; its core = the support level
    at removal (the classic core-decomposition peel). Support of a surviving
    r-clique = number of ALIVE s-cliques (all of whose r-subsets are still alive)
    containing it. We model this exactly with the standard hypergraph peel.
    """
    rcliques = list(all_cliques_of_size(adj, n, r))
    if not rcliques:
        return [], {}
    rc_index = {rc: i for i, rc in enumerate(rcliques)}

    # Every s-clique contributes to the support of each of its r-subsets.
    scliques = list(all_cliques_of_size(adj, n, s))
    # For each s-clique, list its r-subsets (as r-clique indices).
    s_to_r = []           # s_to_r[j] = list of r-clique indices in s-clique j
    r_to_s = defaultdict(list)  # r-clique idx -> list of s-clique idx
    for j, sc in enumerate(scliques):
        subs = []
        for combo in itertools.combinations(sc, r):
            idx = rc_index[combo]   # every r-subset of an s-clique is an r-clique
            subs.append(idx)
        s_to_r.append(subs)
        for idx in subs:
            r_to_s[idx].append(j)

    nrc = len(rcliques)
    s_alive = [True] * len(scliques)
    # support[i] = # alive s-cliques whose every r-subset is alive, containing rc i.
    # Equivalently: # alive s-cliques j with i in s_to_r[j]. An s-clique stays
    # "alive" only while ALL its r-subsets are alive; once any r-subset is peeled,
    # the s-clique is gone (it can no longer be a witness for any survivor).
    support = [0] * nrc
    for i in range(nrc):
        support[i] = len(r_to_s[i])

    core = {}
    alive = [True] * nrc
    removed = 0
    cur_level = 0
    # bucket peel
    import heapq
    heap = [(support[i], i) for i in range(nrc)]
    heapq.heapify(heap)
    while removed < nrc:
        sup, i = heapq.heappop(heap)
        if not alive[i] or sup != support[i]:
            continue
        cur_level = max(cur_level, sup)
        core[rcliques[i]] = cur_level
        alive[i] = False
        removed += 1
        # kill every s-clique through i; decrement its OTHER alive r-subsets.
        for j in r_to_s[i]:
            if not s_alive[j]:
                continue
            s_alive[j] = False
            for k in s_to_r[j]:
                if alive[k] and k != i:
                    support[k] -= 1
                    heapq.heappush(heap, (support[k], k))
    return rcliques, core


def s_adjacency(rcliques, core, adj, n, r, s):
    """Edges between r-cliques that co-occur in some s-clique (share an s-witness)."""
    rc_index = {rc: i for i, rc in enumerate(rcliques)}
    edges = set()
    for sc in all_cliques_of_size(adj, n, s):
        subs = [rc_index[c] for c in itertools.combinations(sc, r)]
        for a in range(len(subs)):
            for b in range(a + 1, len(subs)):
                x, y = subs[a], subs[b]
                if x > y:
                    x, y = y, x
                edges.add((x, y))
    return edges


# ------------------------- canonical nucleus profile -------------------------
def nucleus_profile_direct(rcliques, core, edges):
    """Multiset of (k, component_size) over every distinct core level k.

    At level k: induce the subgraph on r-cliques with core >= k; component_size =
    #r-cliques in each connected component (each contributes 1 here, unit weight).
    """
    n = len(rcliques)
    idx_core = [core[rc] for rc in rcliques]
    # Nuclei are defined for k >= 1 only. r-cliques with core 0 belong to no
    # s-clique, hence to no nucleus (the a_Y engine never materializes them as
    # patterns). Exclude level 0 from the profile so both sides share the same
    # universe (the k>=1 nucleus forest).
    levels = sorted({c for c in set(idx_core) if c >= 1}, reverse=True)
    # adjacency lists
    g = defaultdict(list)
    for x, y in edges:
        g[x].append(y)
        g[y].append(x)
    prof = Counter()
    for k in levels:
        active = [i for i in range(n) if idx_core[i] >= k]
        aset = set(active)
        comp = {}
        cid = 0
        for v in active:
            if v in comp:
                continue
            stack = [v]
            comp[v] = cid
            sz = 0
            while stack:
                u = stack.pop()
                sz += 1
                for w in g[u]:
                    if w in aset and w not in comp:
                        comp[w] = cid
                        stack.append(w)
            prof[(k, sz)] += 1
            cid += 1
    return prof


def gt_merge_tree(rcliques, core, edges):
    """Build the ground-truth elder-rule merge tree at r-clique granularity.

    Independent reimplementation (NOT the a_Y mechanics) used only to self-check
    that merge_tree_profile() faithfully reconstructs the (k,size) profile. Each
    node is one merge component; processing r-cliques core-DESC, an r-clique born
    as a singleton, components sharing an s-edge merge under the elder (highest
    birth, tie larger rep). Returns merge-tree node list (same schema as a_Y CSV).
    """
    n = len(rcliques)
    idx_core = [core[rc] for rc in rcliques]
    active_idx = [i for i in range(n) if idx_core[i] >= 1]
    if not active_idx:
        return []
    g = defaultdict(list)
    aset = set(active_idx)
    for x, y in edges:
        if x in aset and y in aset:
            g[x].append(y)
            g[y].append(x)
    order = sorted(active_idx, key=lambda i: (-idx_core[i], i))

    parent = {}
    comp_size = {}
    birth = {}
    state = {}     # dsu-root -> node id
    nodes = []
    activated = set()

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in order:
        k = idx_core[i]
        parent[i] = i
        comp_size[i] = 1
        birth[i] = k
        nd = {"id": len(nodes), "k_birth": k, "k_death": 0, "parent": -1,
              "size_birth": 1, "size_death": 1}
        nodes.append(nd)
        state[i] = nd["id"]
        activated.add(i)
        # collect distinct roots among already-activated s-neighbors + self
        reps = [find(i)]
        seen = {find(i)}
        for w in g[i]:
            if w in activated:
                r2 = find(w)
                if r2 not in seen:
                    seen.add(r2)
                    reps.append(r2)
        if len(reps) >= 2:
            elder = reps[0]
            for r2 in reps:
                if birth[r2] > birth[elder] or (birth[r2] == birth[elder] and r2 > elder):
                    elder = r2
            elder_node = state[elder]
            for r2 in reps:
                if r2 == elder:
                    continue
                cn = state[r2]
                nodes[cn]["k_death"] = k
                nodes[cn]["size_death"] = comp_size[r2]
                nodes[cn]["parent"] = elder_node
            new_root = elder
            for r2 in reps:
                a = find(new_root)
                b = find(r2)
                if a == b:
                    continue
                if comp_size[a] < comp_size[b]:
                    a, b = b, a
                parent[b] = a
                comp_size[a] += comp_size[b]
                new_root = a
            new_root = find(new_root)
            state[new_root] = elder_node
            birth[new_root] = nodes[elder_node]["k_birth"]
    for r2 in set(find(i) for i in active_idx):
        nid = state[r2]
        if nodes[nid]["parent"] == -1:
            nodes[nid]["size_death"] = comp_size[r2]
    return nodes


def merge_tree_profile(nodes):
    """Reconstruct the (k, component_size) profile from a merge-tree node list.

    nodes: list of dict(id,k_birth,k_death,parent,size_birth,size_death).
    size_birth is in r-clique units. At level k a distinct component is rooted at
    each node n that is ALIVE at k and not yet merged: k_death(n) < k <= k_birth(n)
    (a root has k_death == 0; treat any node with parent==-1 as alive down to its
    own birth only at and below... actually we sweep distinct k values present).

    Component size at level k = sum of size_birth over all nodes m in n's subtree
    with k_birth(m) >= k. We compute it by, for each distinct level k, contracting
    every edge whose child k_death >= k (merged at-or-above k) and summing.
    """
    by_id = {nd["id"]: nd for nd in nodes}
    children = defaultdict(list)
    for nd in nodes:
        if nd["parent"] != -1:
            children[nd["parent"]].append(nd["id"])

    levels = sorted({nd["k_birth"] for nd in nodes}, reverse=True)
    prof = Counter()
    for k in levels:
        # A node is "present" (born) at level k iff k_birth >= k.
        # It is the representative of its component at level k iff it is present
        # and it has NOT merged into its parent at-or-above k, i.e. it is a root
        # OR its own k_death < k (it dies, i.e. merges, strictly below k).
        # Component members of a representative = present nodes reachable by going
        # DOWN into children whose k_death >= k (those merged at-or-above k).
        present = {nd["id"] for nd in nodes if nd["k_birth"] >= k}
        # representative iff present and (root or k_death < k)
        reps = [i for i in present
                if (by_id[i]["parent"] == -1 or by_id[i]["k_death"] < k)]
        for rep in reps:
            # gather subtree members merged at-or-above k
            stack = [rep]
            sz = 0
            while stack:
                u = stack.pop()
                if u not in present:
                    continue
                sz += by_id[u]["size_birth"]
                for c in children[u]:
                    # child merged into u at k_death(c); included iff k_death(c) >= k
                    if by_id[c]["k_death"] >= k and c in present:
                        stack.append(c)
            prof[(k, sz)] += 1
    return prof


# ------------------------------ csv parsing ----------------------------------
def parse_hier_csv(path):
    nodes = []
    with open(path) as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            nodes.append({
                "id": int(parts[0]),
                "k_birth": int(round(float(parts[1]))),
                "k_death": int(round(float(parts[2]))),
                "parent": int(parts[3]),
                "size_birth": int(parts[4]),
                "size_death": int(parts[5]),
            })
    return nodes


def parse_core_dist(stdout):
    dist = Counter()
    for line in stdout.splitlines():
        if line.startswith("core="):
            # core=K count=N
            toks = line.replace("core=", "").replace("count=", "").split()
            k = int(round(float(toks[0])))
            c = int(round(float(toks[1])))
            dist[k] += c
    return dist


# ------------------------------ driver ---------------------------------------
def random_graph(n, p, seed):
    rng = random.Random(seed)
    adj = [set() for _ in range(n)]
    edges = []
    for u in range(n):
        for v in range(u + 1, n):
            if rng.random() < p:
                adj[u].add(v)
                adj[v].add(u)
                edges.append((u, v))
    return adj, edges


def write_edges(path, n, edges):
    with open(path, "w") as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")


def run_instance(binary, n, edges, r, s, tmpdir, no_rmerge):
    gpath = os.path.join(tmpdir, "g.edges")
    hpath = os.path.join(tmpdir, "h.csv")
    write_edges(gpath, n, edges)
    if os.path.exists(hpath):
        os.remove(hpath)
    env = dict(os.environ)
    env["PIVOTER_DUMP_HIER"] = hpath
    if no_rmerge:
        env["SCT_NO_RMERGE"] = "1"
    res = subprocess.run([binary, gpath, str(r), str(s)],
                         capture_output=True, text=True, env=env, timeout=60)
    nodes = parse_hier_csv(hpath) if os.path.exists(hpath) else []
    cdist = parse_core_dist(res.stdout)
    return nodes, cdist, res


def main():
    n_inst = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    binary = sys.argv[2] if len(sys.argv) > 2 else BINARY
    if not os.path.exists(binary):
        print(f"binary not found: {binary}", file=sys.stderr)
        sys.exit(2)

    rs_list = [(2, 3), (2, 4), (3, 4)]
    violations = 0
    core_mismatch = 0
    selfcheck_fail = 0
    tested = 0
    rng = random.Random(12345)

    with tempfile.TemporaryDirectory() as tmp:
        i = 0
        while tested < n_inst:
            i += 1
            n = rng.randint(5, 12)
            p = rng.choice([0.35, 0.45, 0.55, 0.65, 0.75])
            r, s = rng.choice(rs_list)
            no_rmerge = rng.random() < 0.5   # exercise BOTH code paths
            adj, edges = random_graph(n, p, seed=rng.randint(0, 1 << 30))
            if len(edges) == 0:
                continue

            rcliques, core = brute_force_nucleus(adj, n, r, s)
            if not rcliques:
                continue  # no r-cliques: nothing to compare
            edges_s = s_adjacency(rcliques, core, adj, n, r, s)
            gt_prof = nucleus_profile_direct(rcliques, core, edges_s)
            if not gt_prof:
                # no s-cliques at all -> all cores 0; a_Y emits core=0; skip
                # (still counts as a trivial instance but adds no signal)
                continue

            # GT core distribution (weighted by count of r-cliques per core).
            # Restrict to k >= 1: a_Y omits core-0 r-cliques (no s-clique support,
            # never materialized as patterns), so compare on the shared universe.
            gt_cdist = Counter()
            for rc in rcliques:
                if core[rc] >= 1:
                    gt_cdist[core[rc]] += 1

            try:
                nodes, cdist, res = run_instance(binary, n, edges, r, s, tmp, no_rmerge)
            except subprocess.TimeoutExpired:
                print(f"[inst {i}] TIMEOUT n={n} r={r} s={s}")
                violations += 1
                tested += 1
                continue

            tested += 1

            # --- self-check: merge_tree_profile reconstruction is faithful ---
            # Build the GT merge tree with an INDEPENDENT elder-rule implementation
            # and confirm its reconstructed profile equals the direct profile. This
            # validates merge_tree_profile() before we trust it on a_Y's CSV.
            gt_nodes = gt_merge_tree(rcliques, core, edges_s)
            gt_recon = merge_tree_profile(gt_nodes)
            if gt_recon != gt_prof:
                selfcheck_fail += 1
                if selfcheck_fail <= 5:
                    print(f"[inst {i}] SELF-CHECK FAIL (reconstruction bug) n={n} r={r} s={s}")
                    print(f"   direct : {dict(gt_prof)}")
                    print(f"   recon  : {dict(gt_recon)}")

            ay_prof = merge_tree_profile(nodes)

            # a_Y may print a core=0 bucket for degenerate direct regions; the
            # nucleus universe is k>=1, so compare on k>=1 only.
            cdist_pos = Counter({k: v for k, v in cdist.items() if k >= 1})
            ok_core = (gt_cdist == cdist_pos)
            ok_prof = (gt_prof == ay_prof)
            if not ok_core:
                core_mismatch += 1
                if core_mismatch <= 5:
                    print(f"[inst {i}] CORE MISMATCH n={n} r={r} s={s} rmerge={'off' if no_rmerge else 'on'}")
                    print(f"   GT  cores: {dict(sorted(gt_cdist.items()))}")
                    print(f"   a_Y cores: {dict(sorted(cdist_pos.items()))}")
            if not ok_prof:
                violations += 1
                if violations <= 8:
                    print(f"[inst {i}] HIER MISMATCH n={n} r={r} s={s} rmerge={'off' if no_rmerge else 'on'} p={p}")
                    only_gt = gt_prof - ay_prof
                    only_ay = ay_prof - gt_prof
                    print(f"   #rcliques={len(rcliques)} edges={edges}")
                    print(f"   only in GT : {dict(only_gt)}")
                    print(f"   only in a_Y: {dict(only_ay)}")

    print()
    print(f"TESTED: {tested} instances  (r,s) in {rs_list}, n in [5,12]")
    print(f"CORE MISMATCHES:      {core_mismatch}")
    print(f"RECON SELF-CHECK FAIL:{selfcheck_fail}")
    print(f"HIERARCHY VIOLATIONS: {violations}")
    if violations == 0 and core_mismatch == 0 and selfcheck_fail == 0:
        print("RESULT: PASS (0 violations, cores agree)")
        sys.exit(0)
    else:
        print("RESULT: FAIL")
        sys.exit(1)


if __name__ == "__main__":
    main()
