#!/usr/bin/env python3
# Auto-generated for 5500842

STUDENT_ID = "5500842"
STUDENT_NAME = "Wentao Jiang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(graph, core_threshold):
        num_nodes = graph.vertex_num
        diversity_counts = [0] * num_nodes

        for center in range(num_nodes):
            neighborhood = graph.adj_list[center]

            # 剪枝：邻居数量不足以形成 k-core
            if len(neighborhood) < core_threshold:
                continue

            diversity_counts[center] = kCoreBaseStructuralDiversity._evaluate_diverse_cores(
                graph, neighborhood, core_threshold
            )

        return diversity_counts

    @staticmethod
    def _evaluate_diverse_cores(graph, neighbors, k):
        sub_nodes = set(neighbors)
        local_graph = {}
        for u in sub_nodes:
            local_graph[u] = set()
            for v in graph.adj_list[u]:
                if v in sub_nodes:
                    local_graph[u].add(v)

        return kCoreBaseStructuralDiversity._count_k_core_clusters(local_graph, k)

    @staticmethod
    def _count_k_core_clusters(subgraph, k):
        if not subgraph:
            return 0

        node_degrees = {u: len(subgraph[u]) for u in subgraph}
        excluded = set()
        q = deque(u for u in node_degrees if node_degrees[u] < k)

        for u in q:
            excluded.add(u)

        while q:
            u = q.popleft()
            for v in subgraph[u]:
                if v not in excluded:
                    node_degrees[v] -= 1
                    if node_degrees[v] < k:
                        excluded.add(v)
                        q.append(v)

        valid_nodes = [u for u in subgraph if u not in excluded]
        if not valid_nodes:
            return 0

        def iterative_dfs(start, seen):
            stack = [start]
            while stack:
                node = stack.pop()
                if node not in seen:
                    seen.add(node)
                    for neighbor in subgraph[node]:
                        if neighbor not in seen and neighbor not in excluded:
                            stack.append(neighbor)

        seen = set()
        components = 0
        for u in valid_nodes:
            if u not in seen:
                iterative_dfs(u, seen)
                components += 1

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
