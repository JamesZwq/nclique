#!/usr/bin/env python3
# Auto-generated for 5218676

STUDENT_ID = "5218676"
STUDENT_NAME = "Josh Quail"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import defaultdict, deque
from heapq import *
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def core_decomposition(G) -> list[int]:
        core = list(map(len, G.adj_list)) # core starts as just the degree
        u_ordered = sorted(range(G.vertex_num), key=core.__getitem__)

        # set up u_pos and bins
        u_pos = {}
        bins = [0]
        max_deg = 0
        for i, u in enumerate(u_ordered):
            u_pos[u] = i
            if core[u] > max_deg:
                bins.extend([i] * (core[u] - max_deg))
                max_deg = core[u]

        # calculate core numbers
        for u in u_ordered:
            for w in G.adj_list[u]:
                if core[w] > core[u]:
                    w_pos = u_pos[w]
                    w_bin = bins[core[w]]
                    if u != w_bin:
                        u_ordered[w_pos] = u_ordered[w_bin]
                        u_ordered[w_bin] = w
                        u_pos[w] = w_bin
                        u_pos[u_ordered[w_bin]] = w_pos
                    w_bin += 1
                    core[w] -= 1

        return core

    @staticmethod
    def count_components(adj_set: dict[int, set[int]], vertex_filter: set[int]) -> int:
        visited = set()
        components = 0

        for n in vertex_filter:
            if n not in visited:
                stack = [n]
                visited.add(n)
                while len(stack) > 0:
                    u = stack.pop()
                    visited.add(u)
                    for v in adj_set[u]:
                        if v in vertex_filter and v not in visited:
                            stack.append(v)
                components += 1

        return components

    @staticmethod
    def process_one(G, N_v: set[int], k: int) -> int:
        # create subgraph nbr_v using adjacency set
        nbr_v = {}
        for u in N_v:
            nbr_v[u] = set(w for w in G.adj_list[u] if w in N_v)

        # add vertices to q with degree < k
        q = deque()
        degree = defaultdict(int)
        for u in N_v:
            for w in nbr_v[u]:
                degree[u] += 1
            if degree[u] < k:
                q.appendleft(u)

        # remove vertices in q, adding neighbours of removed vertex to q if they now also have degree < k
        while len(q) > 0:
            u = q.pop()
            N_v.remove(u) # this gets N_v closer to the set of vertices of K_k(nbr_v)
            for w in nbr_v[u]:
                if w in N_v and degree[w] >= k:
                    degree[w] -= 1
                    if degree[w] < k:
                        q.appendleft(w)

        # after removing all vertices not in the k-cores, N_v has become the set of vertices of K_k(nbr_v)
        return kCoreBaseStructuralDiversity.count_components(nbr_v, N_v)

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

        # core decomposition using bins
        core = kCoreBaseStructuralDiversity.core_decomposition(G)

        # calculate τs
        τs = []
        for v in range(G.vertex_num):
            N_v = set(w for w in G.adj_list[v] if core[w] > k)
            τ = kCoreBaseStructuralDiversity.process_one(G, N_v, k)
            τs.append(τ)
        return τs

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
