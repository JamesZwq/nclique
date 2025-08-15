#!/usr/bin/env python3
# Auto-generated for 5508494

STUDENT_ID = "5508494"
STUDENT_NAME = "Hrishita Sunil Patil"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity:
    @staticmethod
    def process(G, k):
        n = G.vertex_num  # Total number of vertices in the graph

        # Convert adjacency list to sets for O(1) neighbor lookups
        adj = [set(neigh) for neigh in G.adj_list]

        # Initialize the τ_k vector with zeros
        tau = [0] * n

        # Iterate over each vertex v in the graph
        for v in range(n):
            neighbors = adj[v]  # Get neighbors of vertex v

            # If v has fewer than k neighbors, it can't form a k-core
            if len(neighbors) < k:
                continue

            # Build the induced subgraph on v's neighbors
            subgraph = {}  # Stores induced adjacency list of each neighbor
            deg = {}       # Degree of each neighbor in the induced subgraph

            for u in neighbors:
                common = adj[u] & neighbors  # Induced neighbors (intersection)
                subgraph[u] = common         # Save induced neighbor list
                deg[u] = len(common)         # Degree in the induced subgraph

            # Perform k-core peeling (remove nodes with degree < k)
            q = deque([u for u in deg if deg[u] < k])  # Start with low-degree nodes

            while q:
                u = q.popleft()
                for w in subgraph[u]:
                    if deg[w] >= k:
                        deg[w] -= 1  # Update degree
                        if deg[w] == k - 1:
                            q.append(w)  # Re-check this node
                deg[u] = -1  # Mark u as removed from subgraph

            # Extract remaining nodes as k-core
            k_core_nodes = {u for u in deg if deg[u] >= k}
            if not k_core_nodes:
                continue  # No valid k-core component

            # Count connected components in the k-core
            visited = set()
            count = 0  # Number of connected components in k-core

            for u in k_core_nodes:
                if u not in visited:
                    # Start BFS from this unvisited node
                    count += 1
                    dq = deque([u])
                    visited.add(u)

                    while dq:
                        curr = dq.popleft()
                        for nb in subgraph[curr]:
                            if nb in k_core_nodes and nb not in visited:
                                visited.add(nb)
                                dq.append(nb)

            # Save τ_k(v): number of k-core components in v's neighbor-induced subgraph
            tau[v] = count

        # Return τ_k vector for all vertices
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
