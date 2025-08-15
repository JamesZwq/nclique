#!/usr/bin/env python3
# Auto-generated for 5570877

STUDENT_ID = "5570877"
STUDENT_NAME = "Ming Gao"

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
            if len(neighbors) < k:
                continue

            subgraph = kCoreBaseStructuralDiversity._build_induced_subgraph(G, neighbors)

            kcore_nodes = kCoreBaseStructuralDiversity._extract_k_core_nodes(subgraph, neighbors, k)
            if len(kcore_nodes) < k:
                continue

            sub_kcore_subgraph = {u: subgraph[u] & kcore_nodes for u in kcore_nodes}
            components = kCoreBaseStructuralDiversity._get_connected_components(sub_kcore_subgraph)

            sd[v] = len(components)

        return sd

    @staticmethod
    def _build_induced_subgraph(G, nodes):
        node_set = set(nodes)
        subgraph = {u: set() for u in node_set}
        for u in node_set:
            for v in G.adj_list[u]:
                if v in node_set:
                    subgraph[u].add(v)
        return subgraph

    @staticmethod
    def _extract_k_core_nodes(subgraph, nodes, k):
        subgraph_copy = {u: set(subgraph[u]) & set(nodes) for u in nodes}
        degrees = {u: len(subgraph_copy[u]) for u in nodes}
        queue = deque([u for u in nodes if degrees[u] < k])

        while queue:
            u = queue.popleft()
            for v in subgraph_copy[u]:
                if degrees[v] > 0:
                    degrees[v] -= 1
                    subgraph_copy[v].discard(u)
                    if degrees[v] == k - 1:
                        queue.append(v)
            degrees[u] = -1

        return {u for u in nodes if degrees[u] >= k}

    @staticmethod
    def _get_connected_components(subgraph):
        visited = set()
        components = []

        for node in subgraph:
            if node not in visited:
                comp = []
                queue = deque([node])
                visited.add(node)
                while queue:
                    u = queue.popleft()
                    comp.append(u)
                    for v in subgraph[u]:
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)
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
