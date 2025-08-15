#!/usr/bin/env python3
# Auto-generated for 5190542

STUDENT_ID = "5190542"
STUDENT_NAME = "Courtney Zhang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                continue

            neighbor_set = set(neighbors)
            subgraph = {}
            for u in neighbor_set:
                subgraph[u] = []
                for w in G.adj_list[u]:
                    if w in neighbor_set:
                        subgraph[u].append(w)

            sd[v] = kCoreBaseStructuralDiversity._compute_k_core_count(subgraph, k)

        return sd

    @staticmethod
    def _compute_k_core_count(subgraph, k):
        if not subgraph:
            return 0

        degrees = {u: len(adj) for u, adj in subgraph.items()}

        to_remove = set()
        queue = deque()

        for vertex in subgraph:
            if degrees[vertex] < k:
                to_remove.add(vertex)
                queue.append(vertex)

        while queue:
            current = queue.popleft()
            for neighbor in subgraph[current]:
                if neighbor not in to_remove:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        to_remove.add(neighbor)
                        queue.append(neighbor)

        core_vertices = [v for v in subgraph if v not in to_remove]

        if not core_vertices:
            return 0

        return kCoreBaseStructuralDiversity._count_components(subgraph, core_vertices)

    @staticmethod
    def _count_components(graph, valid_vertices):
        valid_set = set(valid_vertices)
        visited = set()
        components = 0

        for start in valid_vertices:
            if start in visited:
                continue

            components += 1
            bfs_queue = deque([start])
            visited.add(start)

            while bfs_queue:
                node = bfs_queue.popleft()
                for adj in graph[node]:
                    if adj in valid_set and adj not in visited:
                        visited.add(adj)
                        bfs_queue.append(adj)

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
