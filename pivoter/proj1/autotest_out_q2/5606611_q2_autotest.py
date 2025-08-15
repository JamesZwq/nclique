#!/usr/bin/env python3
# Auto-generated for 5606611

STUDENT_ID = "5606611"
STUDENT_NAME = "Tengyao Liu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity:
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

            # Build induced subgraph for neighbors of v
            induced_nodes = set(neighbors)
            induced_subgraph = {}
            for u in induced_nodes:
                induced_subgraph[u] = [x for x in G.adj_list[u] if x in induced_nodes]

            # Remove nodes that can not stay in k-core
            remaining_nodes = kCoreBaseStructuralDiversity.prune_k_core(induced_subgraph, k)

            if not remaining_nodes:
                continue

            # Count components in the remaining subgraph
            count = kCoreBaseStructuralDiversity.count_components(induced_subgraph, remaining_nodes)
            sd[v] = count

        return sd

    @staticmethod
    def prune_k_core(graph, k):
        # Remove nodes with degree < k until stable
        degrees = {u: len(graph[u]) for u in graph}
        removed_nodes = set()
        queue = deque([u for u, d in degrees.items() if d < k])

        for u in queue:
            removed_nodes.add(u)

        while queue:
            node = queue.popleft()
            for neighbor in graph[node]:
                if neighbor not in removed_nodes:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        removed_nodes.add(neighbor)
                        queue.append(neighbor)
        # Returns the set of remaining nodes
        return {u for u in graph if u not in removed_nodes}

    @staticmethod
    def count_components(graph, nodes):
        # Count connected components in the remaining nodes.
        visited = set()
        component_count = 0

        for start in nodes:
            if start in visited:
                continue
            component_count += 1
            stack = [start]
            visited.add(start)

            while stack:
                current = stack.pop()
                for neighbor in graph[current]:
                    if neighbor in nodes and neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)

        return component_count

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
