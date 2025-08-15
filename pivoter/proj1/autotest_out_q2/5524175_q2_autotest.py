#!/usr/bin/env python3
# Auto-generated for 5524175

STUDENT_ID = "5524175"
STUDENT_NAME = "Jingwen Wu"

# ======= 学生代码 =======
from collections import deque
from typing import List, Set


class kCoreBaseStructuralDiversity:
    """
    Compute τ_k(v) for every vertex v in an undirected, unweighted graph:

        τ_k(v) = number of connected components
                 in the k-core(s) of the neighbour-induced subgraph G[N(v)]
    """

    # ------------------------------------------------------------------ #
    # 0. Build an adjacency list from different possible graph APIs      #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _build_adj(G) -> List[List[int]]:
        """Return a list-of-lists adjacency representation."""

        # Already a plain adjacency list
        if isinstance(G, list):
            return [list(neis) for neis in G]

        # Typical course object: G.adj_list_out
        if hasattr(G, "adj_list_out"):
            out = G.adj_list_out
            if isinstance(out, list):                # list indexed by vertex id
                return [list(neis) for neis in out]
            # dict  {v : iterable(neis)}
            return [list(out.get(v, ())) for v in range(G.vertex_num)]

        # Alternative name: G.adj_list
        if hasattr(G, "adj_list"):
            al = G.adj_list
            if isinstance(al, list):
                return [list(neis) for neis in al]
            return [list(al.get(v, ())) for v in range(G.vertex_num)]

        # Generic edges() or edges attribute
        if hasattr(G, "edges"):
            n = G.vertex_num
            adj = [[] for _ in range(n)]
            edges = G.edges() if callable(G.edges) else G.edges
            for u, v in edges:
                adj[u].append(v)
                adj[v].append(u)
            return adj

        # neighbours(v) method
        if hasattr(G, "neighbors"):
            return [list(G.neighbors(v)) for v in range(G.vertex_num)]

        raise TypeError(
            "Unsupported graph type: need adj_list_out / adj_list / "
            "edges() / neighbors(v) or a raw adjacency list."
        )

    # ------------------------------------------------------------------ #
    # 1. Linear-time core-number decomposition (Batagelj–Zaversnik)      #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _core_numbers(adj: List[List[int]]) -> List[int]:
        """
        Return an array core[v] = coreness of v.
        Runs in O(|E|).
        """
        n = len(adj)
        deg = [len(nei) for nei in adj]
        if n == 0:
            return []

        # bucket-sort style ordering by current degree
        maxd = max(deg)
        count = [0] * (maxd + 1)
        for d in deg:
            count[d] += 1

        start = [0] * (maxd + 1)
        ptr = 0
        for d in range(maxd + 1):
            start[d] = ptr
            ptr += count[d]

        vert = [0] * n        # position  -> vertex id
        pos  = [0] * n        # vertex id -> position
        next_free = start[:]

        for v, d in enumerate(deg):
            p = next_free[d]
            vert[p] = v
            pos[v]  = p
            next_free[d] += 1

        # Core peeling
        for i in range(n):
            v = vert[i]
            for u in adj[v]:
                if deg[u] > deg[v]:                # only peel "above" v
                    du = deg[u]
                    pu = pos[u]
                    pw = start[du]
                    w  = vert[pw]
                    if u != w:
                        # swap u with bucket head w
                        vert[pu], vert[pw] = vert[pw], vert[pu]
                        pos[u],  pos[w]   = pw,      pu
                    start[du] += 1
                    deg[u]   -= 1
        return deg                                  # final coreness array

    # ------------------------------------------------------------------ #
    # 2. On a given vertex set S, peel its k-core and count components   #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _k_core_component_count(S: List[int],
                                adj: List[List[int]],
                                k: int) -> int:
        """Return #connected components of the k-core inside S."""
        if len(S) < k:
            return 0

        # initial degree within S
        d = {v: 0 for v in S}
        for v in S:
            d[v] = sum(1 for u in adj[v] if u in d)

        # peel vertices with degree < k
        q = deque(v for v in S if d[v] < k)
        alive: Set[int] = set(S)
        while q:
            v = q.popleft()
            if v not in alive:
                continue
            alive.remove(v)
            for u in adj[v]:
                if u in alive:
                    d[u] -= 1
                    if d[u] == k - 1:
                        q.append(u)

        # count connected components with DFS
        comp = 0
        stack = []
        while alive:
            comp += 1
            v = alive.pop()
            stack.append(v)
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y in alive:
                        alive.remove(y)
                        stack.append(y)
        return comp

    # ------------------------------------------------------------------ #
    # 3. Public entry point                                              #
    # ------------------------------------------------------------------ #
    @staticmethod
    def process(G, k: int) -> List[int]:
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph (any compatible object, see _build_adj)
        k : non-negative integer

        Returns
        -------
        tau : list[int] where tau[v] == τ_k(v)
        """
        n   = G.vertex_num
        adj = kCoreBaseStructuralDiversity._build_adj(G)

        # k == 0 : τ_0(v) == |N(v)|   (every singleton vertex is a 0-core)
        if k <= 0:
            return [len(adj[v]) for v in range(n)]

        # global core numbers, used for pruning
        core = kCoreBaseStructuralDiversity._core_numbers(adj)
        qualifies = [core[v] >= k for v in range(n)]

        tau = [0] * n
        for v in range(n):
            if len(adj[v]) < k:                       # impossible to host k-core
                continue

            # filter neighbours whose global core < k
            nbr = [u for u in adj[v] if qualifies[u]]
            if len(nbr) < k:
                continue

            tau[v] = kCoreBaseStructuralDiversity._k_core_component_count(
                nbr, adj, k
            )

        # cast to plain Python int (avoid numpy scalars)
        return [int(x) for x in tau]

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
