#!/usr/bin/env python3
# Auto-generated for 5511252

STUDENT_ID = "5511252"
STUDENT_NAME = "Chang Liu"

# ======= 学生代码 =======
from collections import deque
from typing import List


class kCoreBaseStructuralDiversity(object):
    """
    Compute τ_k(v) – the k-core-based structural diversity of every vertex
    for an undirected, unweighted graph.

    Key idea :
        1. Global core decomposition (O(|V|+|E|)).
        2. Keep only vertices whose global core ≥ k+1
           – they form the (k+1)-core sub-graph.
        3. Identify each connected component in that sub-graph.
        4. For every component, batch-process the vertices that can
           possibly be affected by it (the component itself + its
           immediate neighbours), and update τ_k.
    """

    # ============================================================= #
    #                           entry                           #
    # ============================================================= #
    @staticmethod
    def process(G, k: int) -> List[int]:
        """
        
        G : UndirectedUnweightedGraph
            
            Required fields:
                • vertex_num : int
                • adj_list   : List[List[int]]  or  Dict[int, List[int]]

        k : int
            The “k” in k-core (k ≥ 0).

        Returns
        -------
        List[int]
            τ_k(v) for every vertex v, indexed by vertex id 0 … |V|-1.
        """
        k = max(k, 0)

        # ---------- basic info ----------
        n = G.vertex_num                       # guaranteed to include isolated vertices
        raw_adj = G.adj_list                   # list or dict – unify to a list below

        # build an adjacency list 'adj' of length n
        if isinstance(raw_adj, list):
            adj = raw_adj + [[] for _ in range(n - len(raw_adj))]
        else:
            # dict case
            adj = [raw_adj.get(i, []) for i in range(n)]

        # ---------- 1. global core numbers ----------
        core = _compute_core_numbers(adj, n)   # O(|V|+|E|) linear peeling

        # ---------- 2. (k+1)-core sub-graph ----------
        k1_core_vertices = [v for v in range(n) if core[v] >= k + 1]
        if not k1_core_vertices:
            return [0] * n  # no structure large enough

        k1_components = _find_components(set(k1_core_vertices), adj)

        # ---------- 3. process every component ----------
        sd = [0] * n                           # τ_k result array
        for comp in k1_components:
            _process_component_batch(comp, adj, core, k, sd)

        return sd


# ==================================================================== #
#                           helper functions                           #
# ==================================================================== #
def _compute_core_numbers(adj: List[List[int]], n: int) -> List[int]:
    """
    Linear-time global core decomposition for an undirected graph.

    Returns List[int]
        core[v] = max k such that v belongs to the k-core of G.
    """
    deg = [len(adj[v]) for v in range(n)]
    if n == 0:
        return []

    # counting-sort buckets by degree
    max_deg = max(deg)
    cnt = [0] * (max_deg + 1)
    for d in deg:
        cnt[d] += 1

    start = [0] * (max_deg + 1)   # starting index of each bucket
    s = 0
    for d in range(max_deg + 1):
        start[d], s = s, s + cnt[d]

    vert = [0] * n     # vertices ordered by non-decreasing degree
    pos = [0] * n      # vert index of each vertex
    next_slot = start[:]  # next free position in each bucket
    for v, d in enumerate(deg):
        idx = next_slot[d]
        vert[idx] = v
        pos[v] = idx
        next_slot[d] += 1

    core = [0] * n
    for i in range(n):
        v = vert[i]
        core[v] = deg[v]
        for u in adj[v]:
            if deg[u] > deg[v]:
                # move u one bucket left (degree –1)
                du = deg[u]
                pu = pos[u]
                pw = start[du]
                w = vert[pw]

                if u != w:        # swap u and the first vertex in its bucket
                    vert[pu], vert[pw] = w, u
                    pos[u], pos[w] = pw, pu

                start[du] += 1
                deg[u] -= 1
    return core


def _find_components(vertices: set, adj: List[List[int]]) -> List[List[int]]:
    """
    Breadth-first search to enumerate connected components
    inside a given vertex set.

    Returns
    -------
    List[List[int]]
        A list of components, each component is a list of vertices.
    """
    visited, comps = set(), []
    for v in vertices:
        if v in visited:
            continue
        comp, dq = [], deque([v])
        visited.add(v)
        while dq:
            u = dq.popleft()
            comp.append(u)
            for w in adj[u]:
                if w in vertices and w not in visited:
                    visited.add(w)
                    dq.append(w)
        comps.append(comp)
    return comps


def _process_component_batch(
    k1_comp: List[int],
    adj: List[List[int]],
    core: List[int],
    k: int,
    sd: List[int],
):
    """
    Given one (k+1)-core connected component, update τ_k for
    all vertices that could possibly be affected by it:
        – vertices inside the component, and
        – their immediate neighbours.

    For each such vertex v:
        – keep only neighbours whose global core ≥ k
          (they are the only ones that can survive a k-core).
        – compute how many k-cores exist in that induced sub-graph.
    """
    comp_set = set(k1_comp)

    # gather the “relevant” vertices to inspect
    relevant = set(k1_comp)
    for v in k1_comp:
        relevant.update(adj[v])

    for v in relevant:
        # effective neighbours = neighbours with core ≥ k
        eff_nbrs = [u for u in adj[v] if core[u] >= k]
        if len(eff_nbrs) < k + 1:
            continue  # impossible to form a k-core

        k_core_cnt = _compute_k_core_count(eff_nbrs, adj, k)
        sd[v] = max(sd[v], k_core_cnt)  # take the max if v is updated multiple times


def _compute_k_core_count(nodes: List[int], adj: List[List[int]], k: int) -> int:
    """
    compute the number of connected k-cores
    in the sub-graph induced by 'nodes'.

    Returns
    -------
    int
        the # of k-cores (i.e., connected components after peeling).
    """
    if len(nodes) < k + 1:
        return 0
    S = set(nodes)

    # initialise degrees inside S
    deg = {u: 0 for u in S}
    for u in S:
        for w in adj[u]:
            if w in S:
                deg[u] += 1

    # peel vertices with degree < k
    q = deque([u for u, d in deg.items() if d < k])
    removed = set()
    while q:
        u = q.popleft()
        if u in removed:
            continue
        removed.add(u)
        for w in adj[u]:
            if w in S and w not in removed:
                deg[w] -= 1
                if deg[w] < k:
                    q.append(w)

    remaining = S - removed
    if not remaining:
        return 0
    return _count_components_set(remaining, adj)


def _count_components_set(vertices: set, adj: List[List[int]]) -> int:
    """
    Depth-first search to count connected components
    inside a vertex *set* (edges restricted to the set).
    """
    visited, comp_cnt = set(), 0
    for v in vertices:
        if v in visited:
            continue
        comp_cnt += 1
        stack = [v]
        visited.add(v)
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w in vertices and w not in visited:
                    visited.add(w)
                    stack.append(w)
    return comp_cnt

# ======= 测试框架 =======

import glob, os, re, time

class UndirectedUnweightedGraph:
    def __init__(self, edge_list):
        info = True
        for u, v in edge_list:
            if info:
                info = False
                self.vertex_num, self.edge_num = u, v
                self.adj_list = [list() for _ in range(self.vertex_num)]
            else:
                self.adj_list[u].append(v)
                self.adj_list[v].append(u)

# 数据集根目录，请按需修改 BASE_DIR
BASE_DIR = "./COMP9312-25T2-Project"

def load_graph(path):
    edges = []
    with open(path) as f:
        for line in f:
            u, v = map(int, line.split())
            edges.append([u, v])
    return UndirectedUnweightedGraph(edges)

def load_expected(path, k_in_name):
    nums = list(map(int, open(path).read().strip().split()))
    return nums[1:] if nums and nums[0] == k_in_name else nums

def run_all_tests():
    start_time = time.time()
    graph_pat  = re.compile(r"data_(\d+)\.graph\.txt")
    answer_pat = re.compile(r"ans_(\d+)_(\d+)\.txt")

    total_expected = 0
    total_mismatch = 0

    for gfile in sorted(glob.glob(os.path.join(BASE_DIR, "data_*.graph.txt"))):
        G = load_graph(gfile)
        size = graph_pat.match(os.path.basename(gfile)).group(1)
        for afile in sorted(glob.glob(os.path.join(BASE_DIR, f"ans_{size}_*.txt"))):
            k = int(answer_pat.match(os.path.basename(afile)).group(2))
            expected = load_expected(afile, k)

            start = time.time()
            tau = kCoreBaseStructuralDiversity.process(G, k)
            # 不打印中间信息
            tau.sort()
            counts = [0] * (tau[-1]+1) if tau else [0]
            for t in tau:
                counts[t] += 1

            # 统计
            total_expected += sum(expected)
            total_mismatch += sum(abs(c - e) for c, e in zip(counts, expected))

    total_time = time.time() - start_time
    # 计算正确率和分数
    total_correct = total_expected - total_mismatch
    correct_rate = total_correct / total_expected if total_expected else 0
    score = correct_rate * 6
    # 输出：zid, 姓名, 正确率, 分数, 总时长
    print(f"{STUDENT_ID}\t{STUDENT_NAME}\t{correct_rate:.2%}\t{score:.2f}\t{total_time:.2f}s")

if __name__ == '__main__':
    run_all_tests()
