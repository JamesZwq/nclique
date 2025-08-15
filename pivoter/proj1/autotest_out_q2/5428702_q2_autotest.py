#!/usr/bin/env python3
# Auto-generated for 5428702

STUDENT_ID = "5428702"
STUDENT_NAME = "Haoran Song"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if len(neighbors) == 0:
                sd[v] = 0
                continue

            neighbor_graph = kCoreBaseStructuralDiversity._build_neighbor_subgraph(G, neighbors)

            k_cores_count = kCoreBaseStructuralDiversity._count_k_cores(neighbor_graph, k)
            sd[v] = k_cores_count

        return sd

    @staticmethod
    def _build_neighbor_subgraph(G, neighbors):

        neighbor_set = set(neighbors)
        vertex_map = {v: i for i, v in enumerate(neighbors)}
        subgraph = [[] for _ in range(len(neighbors))]

        for i, u in enumerate(neighbors):
            for v in G.adj_list[u]:
                if v in neighbor_set:
                    subgraph[i].append(vertex_map[v])

        return subgraph

    @staticmethod
    def _count_k_cores(graph, k):
        n = len(graph)
        if n == 0:
            return 0

        degree = [len(graph[v]) for v in range(n)]
        queue = deque()
        removed = [False] * n

        for v in range(n):
            if degree[v] < k:
                queue.append(v)
                removed[v] = True

        while queue:
            u = queue.popleft()
            for v in graph[u]:
                if not removed[v]:
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)
                        removed[v] = True

        k_core_vertices = [v for v in range(n) if not removed[v]]

        if not k_core_vertices:
            return 0

        visited = [False] * n
        k_cores_count = 0

        for v in k_core_vertices:
            if not visited[v]:
                kCoreBaseStructuralDiversity._dfs(graph, v, visited, removed)
                k_cores_count += 1

        return k_cores_count

    @staticmethod
    def _dfs(graph, start, visited, removed):
        stack = [start]
        visited[start] = True

        while stack:
            u = stack.pop()
            for v in graph[u]:
                if not visited[v] and not removed[v]:
                    visited[v] = True
                    stack.append(v)

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
