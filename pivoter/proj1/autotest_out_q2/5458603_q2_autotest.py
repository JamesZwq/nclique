#!/usr/bin/env python3
# Auto-generated for 5458603

STUDENT_ID = "5458603"
STUDENT_NAME = "Hao Ding"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        for idx in range(n):
            neighbors = set(G.adj_list[idx])
            if not neighbors:
                continue

            subgraph = kCoreBaseStructuralDiversity._build_subgraph(G.adj_list, neighbors)
            sd[idx] = len(kCoreBaseStructuralDiversity._extract_k_core_clusters(subgraph, k))

        return sd

    @staticmethod
    def _build_subgraph(adj_map, nodes):

        return {node: [nbr for nbr in adj_map[node] if nbr in nodes] for node in nodes}

    @staticmethod
    def _extract_k_core_clusters(graph, k):

        if not graph:
            return []
        k_core_nodes = kCoreBaseStructuralDiversity._k_core_peeling(graph, k)
        if not k_core_nodes:
            return []
        return kCoreBaseStructuralDiversity._connected_components(graph, k_core_nodes)

    @staticmethod
    def _k_core_peeling(graph, k):

        degree = {v: len(neigh) for v, neigh in graph.items()}
        queue = deque([v for v in graph if degree[v] < k])
        removed = set(queue)

        while queue:
            v = queue.popleft()
            for u in graph[v]:
                if u not in removed:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        removed.add(u)

        return set(graph) - removed

    @staticmethod
    def _connected_components(graph, valid_nodes):

        components = []
        visited = set()

        for node in valid_nodes:
            if node not in visited:
                comp = set()
                queue = deque([node])
                visited.add(node)
                comp.add(node)

                while queue:
                    current = queue.popleft()
                    for neighbor in graph[current]:
                        if neighbor in valid_nodes and neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
                            comp.add(neighbor)

                components.append(comp)

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
