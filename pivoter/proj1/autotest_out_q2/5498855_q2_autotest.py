#!/usr/bin/env python3
# Auto-generated for 5498855

STUDENT_ID = "5498855"
STUDENT_NAME = "Quxi Kuang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from math import inf
################################################################################

class kCoreBaseStructuralDiversity(object):
    """
    Compute τ_k(v) = number of k-core connected components in G[N(v)]

    Public API:
    - process(G, k) -> List[int]
        - G : UndirectedUnweightedGraph  (vertex ids are 0..n-1)
        - k : int (≥ 0)

    Private API:
    - _core_decomposition(G) -> List[int]
        - G : UndirectedUnweightedGraph  (vertex ids are 0..n-1)
    - _local_k_core_and_components(neigh, adj, k) -> int
        - neigh : List[int]                vertices in N(v) with core>=k
        - adj   : full adjacency list
        - k     : current k

    Returns:
    - tau : List[int]
    """

    # Batagelj-Žaveršnik Algorithm: generate the core number for each vertex
    @staticmethod
    def _core_decomposition(G):
        """ Batagelj-Žaveršnik Algorithm: generate the core number for each vertex
        Parameters:
        - G : UndirectedUnweightedGraph  (vertex ids are 0..n-1)

        Returns:
        - core : List[int]
        """
        n = G.vertex_num
        deg = [len(G.adj_list[v]) for v in range(n)] # the degree of each vertex, stored in a list

        if n == 0:
            return []

        max_deg = max(deg)
        bin_cnt = [0] * (max_deg + 1)          # bucket sizes
        for d in deg:
            bin_cnt[d] += 1 # count the number of vertices with each degree

        # prefix sum → starting index of each bucket
        start = 0
        for d in range(max_deg + 1): # for each degree, store the starting index of the bucket
            cnt = bin_cnt[d]
            bin_cnt[d] = start
            start += cnt

        pos = [0] * n # position of v in vert
        vert = [0] * n # inverse array
        for v, d in enumerate(deg): # for each vertex, store the position of the vertex in the bucket
            pos[v] = bin_cnt[d]
            vert[pos[v]] = v
            bin_cnt[d] += 1

        for d in range(max_deg, 0, -1): # shift right 1 to restore starts
            bin_cnt[d] = bin_cnt[d - 1]
        bin_cnt[0] = 0

        core = deg[:] # will be decreased in-place
        for i in range(n): # for each vertex, update the core number
            v = vert[i]
            for u in G.adj_list[v]:
                if core[u] > core[v]:
                    du, pu = core[u], pos[u]
                    pw = bin_cnt[du]
                    w = vert[pw]
                    if u != w: # O(1) swap u ↔ w
                        vert[pu], vert[pw] = vert[pw], vert[pu]
                        pos[u], pos[w] = pos[w], pos[u]
                    bin_cnt[du] += 1
                    core[u] -= 1
        return core

    # Count the number of k-core connected components in G[neigh]
    @staticmethod
    def _local_k_core_and_components(neigh, adj, k):
        """ Count the number of k-core connected components in G[neigh]
        Parameters:
        - neigh : List[int]
        - adj   : full adjacency list
        - k     : current k

        Returns:
        - number of k-core connected components in G[neigh].
        """
        S_set = set(neigh)
        active = {u: True for u in neigh} # mutable membership
        deg = {u: sum(1 for x in adj[u] if x in S_set) for u in neigh}

        # local peeling: remove vertices with degree less than k
        q = deque([u for u in neigh if deg[u] < k])
        while q:
            u = q.popleft()
            if not active[u]:
                continue
            active[u] = False
            for w in adj[u]:
                if w in active and active[w]:
                    deg[w] -= 1
                    if deg[w] < k:
                        q.append(w)

        # count connected components
        comp_cnt = 0
        for u in neigh:
            if not active.get(u, False):
                continue
            comp_cnt += 1
            stack = [u]
            active[u] = False
            while stack: # dfs to count the number of connected components
                x = stack.pop()
                for y in adj[x]:
                    if active.get(y, False):
                        active[y] = False
                        stack.append(y)
        return comp_cnt

    # Public API
    _core_cache = {}     # {id(G) : core_num list}

    @staticmethod
    def process(G, k):
        """Compute τ_k for every vertex.
        Parameters:
        - G : UndirectedUnweightedGraph  (vertex ids are 0..n-1)
        - k : int (≥ 0)

        Returns:
        - tau : List[int]
        """
        n = G.vertex_num
        if n == 0:
            return []

        # lazy core-number cache per graph instance
        gid = id(G)
        if gid not in kCoreBaseStructuralDiversity._core_cache:
            kCoreBaseStructuralDiversity._core_cache[gid] = (kCoreBaseStructuralDiversity._core_decomposition(G))
        core = kCoreBaseStructuralDiversity._core_cache[gid]

        tau = [0] * n

        # per-vertex processing
        for v in range(n):
            # filter neighbors with sufficient global core number
            neigh = [u for u in G.adj_list[v] if core[u] >= k]
            if len(neigh) <= k: # degree too small → no k-core
                continue
            tau[v] = kCoreBaseStructuralDiversity._local_k_core_and_components(neigh, G.adj_list, k)

        return tau

################################################################################
# You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
################################################################################

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
