#!/usr/bin/env python3
# Auto-generated for 5486836

STUDENT_ID = "5486836"
STUDENT_NAME = "Azeem Ahmed Khan Azeem Ahmed Khan"

# ======= 学生代码 =======
from collections import defaultdict

class UndirectedUnweightedGraph:


    def __init__(self, edges):

        if not edges:
            self.vertex_num = 0
            self.adj_list = []
            return

        self.vertex_num = max(max(edge) for edge in edges) + 1
        self.adj_list = [[] for _ in range(self.vertex_num)]

        # Build adjacency list (undirected)
        for u, v in edges:
            self.adj_list[u].append(v)
            self.adj_list[v].append(u)

        # Remove duplicates and sort
        for i in range(self.vertex_num):
            self.adj_list[i] = sorted(set(self.adj_list[i]))


class kCoreBaseStructuralDiversity:
    """Computes k-core based structural diversity for all vertices."""

    @staticmethod
    def process(G, k):
        """
        Compute τ_k(v) for all vertices v in graph G.

        Args:
            G: UndirectedUnweightedGraph
            k: integer k-core parameter

        Returns:
            List[int]: τ_k(v) for each vertex v
        """
        τ = [0] * G.vertex_num

        for v in range(G.vertex_num):
            neighbors = set(G.adj_list[v])
            k_cores = kCoreBaseStructuralDiversity._find_k_cores(neighbors, G, k)
            τ[v] = len(k_cores)

        return τ

    @staticmethod
    def _find_k_cores(neighbors, G, k):
        """Find all k-core components in neighbor-induced subgraph.
        isku zara dekh ko
        """
        if not neighbors:
            return []

        neighbor_list = list(neighbors)
        n = len(neighbor_list)
        vertex_to_idx = {v: i for i, v in enumerate(neighbor_list)}

        # Build subgraph adjacency list and degrees
        subgraph = [[] for _ in range(n)]
        degrees = [0] * n

        for i, u in enumerate(neighbor_list):
            for v in G.adj_list[u]:
                if v in vertex_to_idx:
                    j = vertex_to_idx[v]
                    subgraph[i].append(j)
                    degrees[i] += 1

        # K-core decomposition: remove vertices with degree < k
        removed = [False] * n
        queue = [i for i in range(n) if degrees[i] < k]

        for i in queue:
            removed[i] = True

        while queue:
            u = queue.pop(0)
            for v in subgraph[u]:
                if not removed[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)
                        removed[v] = True

        # Remaining vertices after pruning
        remaining = [i for i in range(n) if not removed[i]]
        if not remaining:
            return []

        remaining_adj = defaultdict(list)
        for i in remaining:
            for j in subgraph[i]:
                if not removed[j]:
                    remaining_adj[i].append(j)

        # DFS to find connected components
        visited = set()
        components = []

        for start in remaining:
            if start not in visited:
                component = []
                stack = [start]

                while stack:
                    node = stack.pop()
                    if node not in visited:
                        visited.add(node)
                        component.append(neighbor_list[node])
                        stack.extend(n for n in remaining_adj[node] if n not in visited)

                if component:
                    components.append(sorted(component))

        return components

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
