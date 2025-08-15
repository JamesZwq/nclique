#!/usr/bin/env python3
# Auto-generated for 5452016

STUDENT_ID = "5452016"
STUDENT_NAME = "Zixin Zhao"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        adj = G.adj_list
        τ = [0] * n  # Result list to store τ_k(v) for each node

        # Go through each node and handle its neighbor-induced subgraph
        for v in range(n):
            neighbors = adj[v]

            # Not enough neighbors to form a k-core, skip early
            if len(neighbors) < k:
                continue

            # Build the neighbor-induced subgraph #
            neighbor_set = set(neighbors)
            sub_adj = {u: [] for u in neighbor_set}

            for u in neighbor_set:
                for w in adj[u]:
                    if w in neighbor_set:
                        sub_adj[u].append(w)

            # Perform k-core pruning #
            degree = {u: len(sub_adj[u]) for u in sub_adj}
            active = {u: True for u in sub_adj}

            # Initially queue all nodes with degree < k
            q = deque([u for u in sub_adj if degree[u] < k])

            while q:
                u = q.popleft()
                if not active[u]:
                    continue
                active[u] = False  # Mark node as removed

                for v_ in sub_adj[u]:
                    if active[v_]:
                        degree[v_] -= 1
                        if degree[v_] == k - 1:
                            q.append(v_)

            # Build the remaining k-core subgraph #
            kcore_adj = {}
            for u in sub_adj:
                if active[u]:
                    kcore_adj[u] = [v_ for v_ in sub_adj[u] if active[v_]]

            # Count the number of connected components in k-core #
            visited = set()
            count = 0

            for node in kcore_adj:
                if node not in visited:
                    kCoreBaseStructuralDiversity.dfs(node, kcore_adj, visited)
                    count += 1

            τ[v] = count  # Final τ_k(v) for node v

        return τ

    # DFS to explore connected components in k-core subgraph
    @staticmethod
    def dfs(start, graph, visited):
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            for neighbor in graph.get(node, []):
                if neighbor not in visited:
                    stack.append(neighbor)

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
