#!/usr/bin/env python3
# Auto-generated for 5529590

STUDENT_ID = "5529590"
STUDENT_NAME = "Xinxin Ji"

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
        result = [0] * G.vertex_num

        if G.vertex_num == 0:
            return result

        #  5000 only use global core if graph is large
        use_global_core = G.vertex_num > 5000
        global_core = None
        if use_global_core:
            global_core = kCoreBaseStructuralDiversity._compute_global_core(G)

        for v in range(G.vertex_num):
            if use_global_core:
                neighbors = {u for u in G.adj_list[v] if global_core[u] >= k}
            else:
                neighbors = set(G.adj_list[v])

            if not neighbors:
                result[v] = 0
                continue

            # Step 1: Build neighbor-induced subgraph
            sub_adj = {}
            for u in neighbors:
                sub_adj[u] = []
                for w in G.adj_list[u]:
                    if w in neighbors:
                        sub_adj[u].append(w)

            # Step 2: Compute k-core in the subgraph using optimized method
            k_core_nodes = kCoreBaseStructuralDiversity._compute_k_core(sub_adj, k)
            if not k_core_nodes:
                result[v] = 0
                continue

            # Step 3: Count connected components
            components = kCoreBaseStructuralDiversity._find_connected_components(sub_adj, k_core_nodes)
            result[v] = len(components)

        return result

    @staticmethod
    def _compute_global_core(G):
        degree = [len(G.adj_list[v]) for v in range(G.vertex_num)]
        core = degree[:]
        changed = True

        while changed:
            changed = False
            for v in range(G.vertex_num):
                if core[v] > degree[v]:
                    core[v] = degree[v]
                    for u in G.adj_list[v]:
                        if core[u] > core[v]:
                            degree[u] -= 1
                            changed = True

        return core

    @staticmethod
    def _compute_k_core(graph, k):
        degree = {v: len(graph[v]) for v in graph}
        queue = deque([v for v in graph if degree[v] < k])
        removed = set(queue)

        while queue:
            u = queue.popleft()
            for neighbor in graph[u]:
                if neighbor not in removed:
                    degree[neighbor] -= 1
                    if degree[neighbor] < k:
                        queue.append(neighbor)
                        removed.add(neighbor)

        return set(graph.keys()) - removed

    @staticmethod
    def _find_connected_components(graph, nodes):
        visited = set()
        components = []

        for v in nodes:
            if v in visited:
                continue

            component = []
            queue = deque([v])
            visited.add(v)

            while queue:
                u = queue.popleft()
                component.append(u)
                for neighbor in graph[u]:
                    if neighbor in nodes and neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)

            components.append(component)

        return components



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
