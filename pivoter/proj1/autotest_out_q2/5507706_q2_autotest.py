#!/usr/bin/env python3
# Auto-generated for 5507706

STUDENT_ID = "5507706"
STUDENT_NAME = "Ariel Sun"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        total_nodes = G.vertex_num
        diversity = [0] * total_nodes

        for node in range(total_nodes):
            neighbour_list = G.adj_list[node]
            if not neighbour_list:
                continue

            # Phase 1: Build subgraph of the node's neighbours
            subgraph = kCoreBaseStructuralDiversity.create_neighbour_subgraph(G, neighbour_list)

            # Phase 2: Prune nodes to get k-core
            remaining_nodes = kCoreBaseStructuralDiversity.prune_k_core(subgraph, k)

            # Phase 3: Count connected components in the pruned subgraph
            if remaining_nodes:
                component_count = kCoreBaseStructuralDiversity.count_components(subgraph, remaining_nodes)
                diversity[node] = component_count

        return diversity

    @staticmethod
    def create_neighbour_subgraph(G, neighbours):
        node_set = set(neighbours)
        subgraph = {}
        for node in node_set:
            connections = []
            for adj in G.adj_list[node]:
                if adj in node_set:
                    connections.append(adj)
            subgraph[node] = connections
        return subgraph

    @staticmethod
    def prune_k_core(graph, k):
        degrees = {u: len(adj) for u, adj in graph.items()}
        queue = deque([u for u in graph if degrees[u] < k])
        removed = set(queue)

        while queue:
            current = queue.popleft()
            for neighbor in graph[current]:
                if neighbor not in removed:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        removed.add(neighbor)
                        queue.append(neighbor)

        # Return the set of nodes that survived the pruning
        return [u for u in graph if u not in removed]

    @staticmethod
    def count_components(graph, nodes):
        visited = set()
        count = 0

        for start in nodes:
            if start in visited:
                continue
            count += 1
            stack = [start]
            visited.add(start)
            while stack:
                current = stack.pop()
                for neighbor in graph[current]:
                    if neighbor not in visited and neighbor in nodes:
                        visited.add(neighbor)
                        stack.append(neighbor)
        return count

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
