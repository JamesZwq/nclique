#!/usr/bin/env python3
# Auto-generated for 5470992

STUDENT_ID = "5470992"
STUDENT_NAME = "Yiqiao Cheng"

# ======= 学生代码 =======
from collections import deque, defaultdict

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        # Helper function to build the neighbor-induced subgraph of a vertex
        def build_induced_subgraph(neighbors):
            subgraph = defaultdict(set)
            neighbor_set = set(neighbors)
            for u in neighbors:
                for v in G.adj_list[u]:
                    if v in neighbor_set:
                        subgraph[u].add(v)
                        subgraph[v].add(u)
            return subgraph

        # Performs k-core decomposition on a given subgraph
        def k_core_decomposition(subgraph, k):
            # Initialize degrees of all nodes in the subgraph
            degrees = {node: len(neigh) for node, neigh in subgraph.items()}
            queue = deque([node for node, deg in degrees.items() if deg < k])

            # Iteratively remove nodes whose degree is less than k
            while queue:
                node = queue.popleft()
                if node not in degrees:
                    continue
                for neighbor in subgraph[node]:
                    if neighbor in degrees:
                        degrees[neighbor] -= 1
                        if degrees[neighbor] < k:
                            queue.append(neighbor)
                del degrees[node]

            # Return remaining nodes that belong to the k-core
            return set(degrees.keys())

        # Count connected components in the k-core subgraph
        def count_connected_components(nodes_in_kcore, subgraph):
            visited = set()
            count = 0
            for node in nodes_in_kcore:
                if node not in visited:
                    count += 1
                    queue = deque([node])
                    visited.add(node)
                    while queue:
                        curr = queue.popleft()
                        for neighbor in subgraph[curr]:
                            if neighbor in nodes_in_kcore and neighbor not in visited:
                                visited.add(neighbor)
                                queue.append(neighbor)
            return count

        n = G.vertex_num
        τ = [0] * n  # Result list storing τ_k(v) for each vertex v

        for v in range(n):
            neighbors = set(G.adj_list[v])
            if not neighbors:
                τ[v] = 0
                continue

            # Build the neighbor-induced subgraph for vertex v
            subgraph = build_induced_subgraph(neighbors)

            # Compute the k-core of the neighbor-induced subgraph
            k_core_nodes = k_core_decomposition(subgraph, k)

            # Count how many connected k-cores exist in the subgraph
            τ[v] = count_connected_components(k_core_nodes, subgraph)

        return τ

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
