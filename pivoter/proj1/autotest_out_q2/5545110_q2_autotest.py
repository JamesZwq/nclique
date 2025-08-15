#!/usr/bin/env python3
# Auto-generated for 5545110

STUDENT_ID = "5545110"
STUDENT_NAME = "Jikun Lyu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################
# For each vertex v, we compute τ_k(v), the number of k-core components in the
# subgraph induced by v's neighbors. We use:
#   - global core number pruning,
#   - local k-core peeling (BFS),
#   - union-find for connected component count.
#
# Time Complexity:
#   - Per vertex v: O(deg(v)^2)
#   - Overall: O(v∈V∑ deg(v)^2)
#   - Space: O(N + M)
#
# The final two test cases complete in approximately 1.5 seconds each
#
# Specific algorithm design and full time/space complexity explanation:
# See Q2.pdf
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

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
        n = G.vertex_num
        if n == 0:
            return []

        core = kCoreBaseStructuralDiversity._core_numbers(G.adj_list)
        adj_sets = [set(neigh) for neigh in G.adj_list]
        sd = [0] * n

        for v in range(n):
            if core[v] < k:
                continue

            neighbors = G.adj_list[v]
            if len(neighbors) < k:
                continue

            S = [u for u in neighbors if core[u] >= k]
            if len(S) < k:
                continue

            sd[v] = kCoreBaseStructuralDiversity._count_kcores(S, adj_sets, k)

        return sd

    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################

    @staticmethod
    def _core_numbers(adj_list):
        n = len(adj_list)
        if n == 0:
            return []

        deg = [len(neigh) for neigh in adj_list]
        max_deg = max(deg) if deg else 0

        bin_boundaries = [0] * (max_deg + 2)
        for d in deg:
            bin_boundaries[d] += 1

        start = 0
        for d in range(max_deg + 1):
            num = bin_boundaries[d]
            bin_boundaries[d] = start
            start += num
        bin_boundaries[max_deg + 1] = n

        pos = [0] * n
        vert = [0] * n
        for v in range(n):
            pos[v] = bin_boundaries[deg[v]]
            vert[pos[v]] = v
            bin_boundaries[deg[v]] += 1

        for d in range(max_deg, 0, -1):
            bin_boundaries[d] = bin_boundaries[d - 1]
        bin_boundaries[0] = 0

        core = deg.copy()
        for i in range(n):
            v = vert[i]
            for u in adj_list[v]:
                if core[u] > core[v]:
                    pu = pos[u]
                    pw = bin_boundaries[core[u]]
                    w = vert[pw]
                    if u != w:
                        pos[u], pos[w] = pw, pu
                        vert[pu], vert[pw] = w, u
                    bin_boundaries[core[u]] += 1
                    core[u] -= 1

        return core

    @staticmethod
    def _count_kcores(nodes, adj_sets, k):
        node_set = set(nodes)
        sub_adj = {u: list(adj_sets[u] & node_set) for u in nodes}
        degrees = {u: len(neigh) for u, neigh in sub_adj.items()}
        q = deque([u for u in degrees if degrees[u] < k])
        removed = set()

        while q:
            u = q.popleft()
            if u in removed:
                continue
            removed.add(u)
            for v in sub_adj[u]:
                if v in degrees and v not in removed:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        q.append(v)

        remaining = node_set - removed
        if not remaining:
            return 0

        parent = {u: u for u in remaining}

        def find(u):
            while parent[u] != u:
                parent[u] = parent[parent[u]]
                u = parent[u]
            return u

        for u in remaining:
            for v in sub_adj[u]:
                if v in remaining and u < v:
                    ru, rv = find(u), find(v)
                    if ru != rv:
                        parent[rv] = ru

        return len(set(find(u) for u in remaining))

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
