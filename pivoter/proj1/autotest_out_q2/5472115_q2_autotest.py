#!/usr/bin/env python3
# Auto-generated for 5472115

STUDENT_ID = "5472115"
STUDENT_NAME = "Asher Yang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute structural diversity (τ) values for all vertices in graph G

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            Graph object with attributes:
            - vertex_num: number of vertices
            - adj_list: adjacency list representation
        k : int
            Minimum degree for k-core

        Returns
        -------
        List[int]
            τ_k(v) for each vertex v (0-indexed)
        """
        τ = [0] * G.vertex_num

        for v in range(G.vertex_num):
            # Step 1: Get neighbors of v
            neighbors = set(G.adj_list[v])
            if not neighbors:
                continue

            # Step 2: Build neighbor-induced subgraph
            neighbor_list = list(neighbors)
            vertex_map = {u: i for i, u in enumerate(neighbor_list)}

            # Create adjacency list for subgraph
            subgraph_adj = [[] for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list):
                for w in G.adj_list[u]:
                    if w in neighbors:
                        subgraph_adj[i].append(vertex_map[w])

            # Step 3: Find k-core in subgraph
            core_vertices = kCoreBaseStructuralDiversity._find_kcore(subgraph_adj, k)

            # Step 4: Count connected components in k-core
            if core_vertices:
                τ[v] = kCoreBaseStructuralDiversity._count_components(subgraph_adj, core_vertices)

        return τ

    @staticmethod
    def _find_kcore(adj_list, k):
        """Find vertices in the k-core using degree pruning"""
        n = len(adj_list)
        degrees = [len(adj) for adj in adj_list]
        active = [True] * n

        # Iteratively prune vertices with degree < k
        changed = True
        while changed:
            changed = False
            for v in range(n):
                # Check if vertex should be pruned
                if active[v] and degrees[v] < k:
                    active[v] = False  # Remove from k-core
                    changed = True
                    # Update degrees of neighbors
                    for u in adj_list[v]:
                        if active[u]:
                            degrees[u] -= 1  # Reduce neighbor's degree

        # Return remaining vertices that satisfy k-core condition
        return {v for v in range(n) if active[v]}

    @staticmethod
    def _count_components(adj_list, vertices):
        """Count connected components in induced subgraph"""
        visited = set()
        components = 0

        for v in vertices:
            if v not in visited:
                # New component found
                components += 1
                # Start DFS from this vertex
                stack = [v]
                visited.add(v)

                while stack:
                    u = stack.pop()
                    # Explore all neighbors in the induced subgraph
                    for w in adj_list[u]:
                        if w in vertices and w not in visited:
                            visited.add(w)
                            stack.append(w)

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
