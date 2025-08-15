#!/usr/bin/env python3
# Auto-generated for 5464397

STUDENT_ID = "5464397"
STUDENT_NAME = "Hongyang Huang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    """
    Calculate τ_k(v): the number of k-cores in the neighbor-induced subgraph of vertex v
    （k-core-based structural diversity）。
    Compatible with two graph interfaces:
        If G provides a getNeighbors(v) method → ​​use it first
        Otherwise the default is to read G.adj_list or G.adj as the adjacency list.
    """

    @staticmethod
    def _neigh_accessor(G):
        """
        Returns a function neigh(u) -> List[int], which is used to get the neighbor list of vertex u.
        If the graph object implements getNeighbors(v), call it directly; otherwise read adj/adj_list.
        """
        if hasattr(G, "getNeighbors") and callable(getattr(G, "getNeighbors")):
            return G.getNeighbors
        else:
            adj = (G.adj_list if hasattr(G, "adj_list") else G.adj)
            return lambda u: adj[u]

    @staticmethod
    def _count_components(neigh, nbr_set, removed):
        """
        Count the number of connected components in nbr_set (after removing removed vertices).
        Use the neighbor access function neigh to avoid relying on specific graph implementation details.
        """
        visited = set()
        comp_cnt = 0
        for u in nbr_set:
            if u in removed or u in visited:
                continue
            comp_cnt += 1
            stack = [u]
            visited.add(u)
            while stack:
                x = stack.pop()
                for y in neigh(x):
                    if y in nbr_set and y not in removed and y not in visited:
                        visited.add(y)
                        stack.append(y)
        return comp_cnt

    @staticmethod
    def process(G, k):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        neigh = kCoreBaseStructuralDiversity._neigh_accessor(G)
        n = G.vertex_num
        tau = [0] * n

        for v in range(n):
            nbrs = neigh(v)
            deg_v = len(nbrs)
            if deg_v == 0 or k > deg_v - 1:
                # No neighbors or k greater than the local maximum possible degree ⇒ no k-core
                continue

            nbr_set = set(nbrs)
            degree = {}
            peel_q = deque()

            # 1. Local Degree + Skinning Initialization
            for u in nbr_set:
                du = sum(1 for w in neigh(u) if w in nbr_set)
                degree[u] = du
                if du < k:
                    peel_q.append(u)

            # 2. k-core peeling
            removed = set(peel_q)
            while peel_q:
                u = peel_q.popleft()
                for w in neigh(u):
                    if w in nbr_set and w not in removed:
                        degree[w] -= 1
                        if degree[w] < k:
                            removed.add(w)
                            peel_q.append(w)

            # 3. Connected component statistics
            tau[v] = kCoreBaseStructuralDiversity._count_components(
                neigh, nbr_set, removed
            )

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
