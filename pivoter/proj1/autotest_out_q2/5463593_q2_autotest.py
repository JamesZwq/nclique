#!/usr/bin/env python3
# Auto-generated for 5463593

STUDENT_ID = "5463593"
STUDENT_NAME = "Wenyuan Zhuang"

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
        # Total number of vertices in the graph
        n = G.vertex_num
        # Initialize structural diversity scores
        sd = [0] * n

        for v in range(n):
            # Get neighbors of vertex v
            neighbors = set(G.adj_list[v])
            # Skip if no neighbors or all neighbors have degree < k
            if not neighbors or all(len(G.adj_list[u]) < k for u in neighbors):
                continue

            # Step 1: Build the neighbor-induced subgraph
            induced_subgraph = kCoreBaseStructuralDiversity.build_neighbor_induced_subgraph(G, neighbors)
            # Step 2: Extract the k-core subgraph from it
            k_core_subgraph = kCoreBaseStructuralDiversity.extract_k_core_subgraph(induced_subgraph, k)
            # Step 3: Count connected components in the k-core subgraph
            sd[v] = kCoreBaseStructuralDiversity.count_connected_components(k_core_subgraph)

        return sd

    @staticmethod
    def build_neighbor_induced_subgraph(G, neighbors):
        """
        Construct an induced subgraph using only the neighbors of a vertex.
        Only retain edges between nodes in the neighbor set.
        """
        adjacency_list = {u: [] for u in neighbors}
        for u in neighbors:
            for w in G.adj_list[u]:
                if w in neighbors:
                    adjacency_list[u].append(w)
        return adjacency_list

    @staticmethod
    def extract_k_core_subgraph(adjacency_list, k):
        """
        Extract the k-core from the given adjacency list.
        Iteratively remove nodes with degree < k.
        """
        degree = {u: len(neighbors) for u, neighbors in adjacency_list.items()}
        to_remove = deque([u for u in degree if degree[u] < k])

        while to_remove:
            u = to_remove.popleft()
            for w in adjacency_list[u]:
                if w in degree:
                    if degree[w] >= k:
                        degree[w] -= 1
                        if degree[w] < k:
                            to_remove.append(w)
            del degree[u]  # Remove node from graph

        # Build adjacency list for the remaining k-core subgraph
        remaining_nodes = set(degree.keys())
        k_core_adj_list = {}
        for u in remaining_nodes:
            k_core_adj_list[u] = []
            for w in adjacency_list[u]:
                if w in remaining_nodes:
                    k_core_adj_list[u].append(w)

        return k_core_adj_list

    @staticmethod
    def count_connected_components(adjacency_list):
        """
        Count the number of connected components in the subgraph
        using depth-first search (DFS).
        """
        visited = set()
        num_components = 0

        for node in adjacency_list:
            if node not in visited:
                num_components += 1
                stack = deque([node])
                visited.add(node)
                while stack:
                    current = stack.pop()
                    for neighbor in adjacency_list[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            stack.append(neighbor)
        return num_components

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
