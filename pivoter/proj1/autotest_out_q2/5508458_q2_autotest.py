#!/usr/bin/env python3
# Auto-generated for 5508458

STUDENT_ID = "5508458"
STUDENT_NAME = "Jinhan Wan"

# ======= 学生代码 =======
from collections import deque, defaultdict

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def compute_k_core_connected_components(subgraph_adj, k):
        """
        Compute the number of connected components in the k-core of the neighbor-induced subgraph.
        """
        degree = {u: len(neighbors) for u, neighbors in subgraph_adj.items()}
        removed = set()

        # Iteratively remove nodes with degree < k
        changed = True
        while changed:
            changed = False
            to_remove = []
            for u in subgraph_adj:
                if u not in removed and degree[u] < k:
                    to_remove.append(u)
                    changed = True
            for u in to_remove:
                removed.add(u)
                for v in subgraph_adj[u]:
                    if v not in removed:
                        degree[v] -= 1

        # Get the remaining nodes after k-core peeling
        remaining_nodes = set(subgraph_adj.keys()) - removed

        # Count the number of connected components in the remaining subgraph
        visited = set()
        def bfs(start):
            queue = deque([start])
            visited.add(start)
            while queue:
                u = queue.popleft()
                for v in subgraph_adj[u]:
                    if v in remaining_nodes and v not in visited:
                        visited.add(v)
                        queue.append(v)

        count = 0
        for node in remaining_nodes:
            if node not in visited:
                bfs(node)
                count += 1

        return count

    @staticmethod
    def process(G, k):
        """
        For each vertex v, compute τ_k(v): the number of connected components in the k-core of its neighbor-induced subgraph.
        """
        n = G.vertex_num
        tau = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                tau[v] = 0
                continue

            # Construct the neighbor-induced subgraph
            sub_nodes = set(neighbors)
            subgraph_adj = defaultdict(list)
            for u in sub_nodes:
                for w in G.adj_list[u]:
                    if w in sub_nodes:
                        subgraph_adj[u].append(w)

            # Compute structural diversity based on k-core
            tau[v] = kCoreBaseStructuralDiversity.compute_k_core_connected_components(subgraph_adj, k)

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
