#!/usr/bin/env python3
# Auto-generated for 5456839

STUDENT_ID = "5456839"
STUDENT_NAME = "Wentao Lu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        tau = [0] * n  # Output list for τ_k(v)

        for v in range(n):
            nbrs = G.adj_list[v]
            if not nbrs:
                continue

            # Build neighbor-induced subgraph of v
            induced_subgraph = {u: set() for u in nbrs}
            for u in nbrs:
                for nei in G.adj_list[u]:
                    if nei in induced_subgraph:
                        induced_subgraph[u].add(nei)

            # Extract k-core from the induced subgraph
            core = kCoreBaseStructuralDiversity._k_core(induced_subgraph, k)

            # Count number of connected components in the k-core
            visited = set()
            components = 0
            for node in core:
                if node not in visited:
                    kCoreBaseStructuralDiversity._dfs(node, core, visited)
                    components += 1

            tau[v] = components  # τ_k(v)

        return tau

    @staticmethod
    def _k_core(graph, k):
        # Iteratively remove nodes with degree < k
        degrees = {u: len(neis) for u, neis in graph.items()}
        changed = True
        while changed:
            changed = False
            for u in list(graph):
                if degrees[u] < k:
                    for v in graph[u]:
                        if v in graph:
                            degrees[v] -= 1
                            graph[v].discard(u)
                    del graph[u]
                    degrees.pop(u)
                    changed = True
        return graph

    @staticmethod
    def _dfs(start, graph, visited):
        # Standard DFS to traverse connected component
        stack = [start]
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                stack.extend(graph[node] - visited)

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
