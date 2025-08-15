#!/usr/bin/env python3
# Auto-generated for 5578805

STUDENT_ID = "5578805"
STUDENT_NAME = "Biao Geng"

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
            neighbors = set(G.adj_list[v])
            if not neighbors:
                sd[v] = 0
                continue

            neighbor_subgraph = kCoreBaseStructuralDiversity.build_neighbor_subgraph(G, neighbors)
            sd[v] = kCoreBaseStructuralDiversity.count_k_cores(neighbor_subgraph, k)

        return sd

    @staticmethod
    def build_neighbor_subgraph(G, neighbors):
        subgraph = {}
        for u in neighbors:
            subgraph[u] = [v for v in G.adj_list[u] if v in neighbors]
        return subgraph


    @staticmethod
    def count_k_cores(adj_dict, k):
        if not adj_dict:
            return 0

        visited = set()
        k_core_count = 0

        for node in adj_dict:
            if node not in visited:
                component = kCoreBaseStructuralDiversity.get_component(adj_dict, node, visited)
                k_cores = kCoreBaseStructuralDiversity.find_k_cores(adj_dict, component, k)
                k_core_count += len(k_cores)

        return k_core_count

    @staticmethod
    def get_component(adj_dict, start, visited):
        component = []
        queue = deque([start])
        visited.add(start)

        while queue:
            u = queue.popleft()
            component.append(u)
            for v in adj_dict.get(u, []):
                if v not in visited:
                    visited.add(v)
                    queue.append(v)
        return component

    @staticmethod
    def find_k_cores(adj_dict, component, k):
        if len(component) < k:
            return []

        comp_set = set(component)
        comp_adj = {u: [v for v in adj_dict[u] if v in comp_set] for u in component}

        k_core_vertices = kCoreBaseStructuralDiversity.k_core_decomposition(comp_adj, k)

        if not k_core_vertices:
            return []

        k_cores = []
        k_core_visited = set()

        for v in k_core_vertices:
            if v not in k_core_visited:
                core_component = []
                queue = deque([v])
                k_core_visited.add(v)

                while queue:
                    u = queue.popleft()
                    core_component.append(u)
                    for neighbor in comp_adj[u]:
                        if neighbor in k_core_vertices and neighbor not in k_core_visited:
                            k_core_visited.add(neighbor)
                            queue.append(neighbor)

                if len(core_component) >= k:
                    k_cores.append(core_component)

        return k_cores

    @staticmethod
    def k_core_decomposition(adj_dict, k):
        degrees = {u: len(adj_dict[u]) for u in adj_dict}
        to_remove = [u for u in degrees if degrees[u] < k]

        while to_remove:
            u = to_remove.pop()
            if u not in degrees:
                continue

            del degrees[u]

            for v in adj_dict[u]:
                if v in degrees:
                    degrees[v] -= 1
                    if degrees[v] < k and v not in to_remove:
                        to_remove.append(v)

        return list(degrees.keys())

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
