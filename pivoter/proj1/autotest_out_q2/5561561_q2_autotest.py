#!/usr/bin/env python3
# Auto-generated for 5561561

STUDENT_ID = "5561561"
STUDENT_NAME = "Lei Shi"

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
        sd = [0] * n  # τ_k(v) for each vertex

        for v in range(n):
            neighbors = G.adj_list[v]
            
            if len(neighbors) < k:
                continue

            subgraph = kCoreBaseStructuralDiversity.get_neighbor_induced_subgraph(G, neighbors)
            core_nodes = kCoreBaseStructuralDiversity.get_nodes(subgraph, k)
            if not core_nodes:
                continue

            core_graph = kCoreBaseStructuralDiversity.build_core_subgraph(subgraph, core_nodes)
            sd[v] = kCoreBaseStructuralDiversity.count_components(core_graph)

        return sd

    @staticmethod
    def get_neighbor_induced_subgraph(G, neighbors):
        neighbors_set = set(neighbors)
        subgraph = dict()
        for u in neighbors:
            subgraph[u] = [w for w in G.adj_list[u] if w in neighbors_set]
        return subgraph

    @staticmethod
    def get_nodes(subgraph, k):
        # remove nodes with degree < k
        degree = {u: len(adj) for u, adj in subgraph.items()}
        queue = deque([u for u, d in degree.items() if d < k])

        while queue:
            u = queue.popleft()
            for v in subgraph[u]:
                if v in degree:
                    degree[v] -= 1
                    if degree[v] == k - 1:
                        queue.append(v)
            del degree[u]

        return set(degree.keys())

    @staticmethod
    def build_core_subgraph(subgraph, core_nodes):
        core_graph = dict()
        for u in core_nodes:
            core_graph[u] = [v for v in subgraph[u] if v in core_nodes]
        return core_graph

    @staticmethod
    def count_components(graph):
        # DFS to count components
        visited = set()
        count = 0

        for node in graph:
            if node not in visited:
                stack = [node]
                while stack:
                    u = stack.pop()
                    if u not in visited:
                        visited.add(u)
                        stack.extend(graph[u])
                count += 1

        return count




    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################

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
