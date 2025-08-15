#!/usr/bin/env python3
# Auto-generated for 5341069

STUDENT_ID = "5341069"
STUDENT_NAME = "Haotian Luo"

# ======= 学生代码 =======
from collections import deque, defaultdict

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        # Build neighbor-induced subgraph (adjacency list)
        def build_induced_subgraph(neighbors):
            sub_adj = defaultdict(set)
            for u in neighbors:
                for v in G.adj_list[u]:
                    if v in neighbors:
                        sub_adj[u].add(v)
            return sub_adj

        # Extract k-core subgraph using degree-based pruning
        def extract_k_core(adj):
            degrees = {node: len(adj[node]) for node in adj}
            queue = deque([node for node in degrees if degrees[node] < k])
            visited = set()

            while queue:
                node = queue.popleft()
                if node in visited:
                    continue
                visited.add(node)

                for neighbor in adj[node]:
                    if neighbor in degrees:
                        degrees[neighbor] -= 1
                        if degrees[neighbor] < k:
                            queue.append(neighbor)

                if node in degrees:
                    del degrees[node]

            # Build new adjacency list for remaining k-core nodes
            k_core_adj = defaultdict(set)
            for u in degrees:
                for v in adj[u]:
                    if v in degrees:
                        k_core_adj[u].add(v)
            return k_core_adj

        # Count the number of connected components using BFS
        def count_connected_components(adj):
            visited = set()
            count = 0
            for node in adj:
                if node not in visited:
                    queue = deque([node])
                    visited.add(node)
                    while queue:
                        u = queue.popleft()
                        for v in adj[u]:
                            if v not in visited:
                                visited.add(v)
                                queue.append(v)
                    count += 1
            return count

        n = G.vertex_num  # Number of nodes
        tau = [0] * n     # Result vector: τ_k(v) for each node

        for v in range(n):
            neighbors = set(G.adj_list[v])
            if not neighbors:
                tau[v] = 0  # No neighbors → diversity is 0
                continue

            sub_adj = build_induced_subgraph(neighbors)     # Step 1: get induced subgraph
            k_core_adj = extract_k_core(sub_adj)            # Step 2: extract k-core
            tau[v] = count_connected_components(k_core_adj) # Step 3: count components

        return tau

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
