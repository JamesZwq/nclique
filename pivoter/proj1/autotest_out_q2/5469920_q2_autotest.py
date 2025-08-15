#!/usr/bin/env python3
# Auto-generated for 5469920

STUDENT_ID = "5469920"
STUDENT_NAME = "Qirui Zhang"

# ======= 学生代码 =======
from collections import deque
from typing import List

class kCoreBaseStructuralDiversity(object):

    @staticmethod
    def _get_adj_list(G):
        # Extract an undirected graph adjacency list from a graph object G
        n = G.vertex_num
        if hasattr(G, "adj_list"):
            return G.adj_list
        if hasattr(G, "adj"):
            return G.adj
        adj = [[] for _ in range(n)]
        if hasattr(G, "edges"):
            for u, v in G.edges:
                adj[u].append(v)
                adj[v].append(u)
        return adj

    @staticmethod
    def process(G, k):
        # Compute the k-core decomposition of a graph G
        # Input: G - a graph object； k - the k-core value
        # Output: τ - a list of integers, where τ[v] is the k-core value of vertex v
        n   = G.vertex_num
        adj = kCoreBaseStructuralDiversity._get_adj_list(G)
        τ   = [0] * n 

        # Boundary 1: k <= 0
        if k <= 0:
            for v in range(n):
                nbrs = adj[v]
                if not nbrs:
                    continue
                seen = set()
                comps = 0
                for u in nbrs:
                    if u in seen:
                        continue
                    comps += 1
                    stack = [u]
                    seen.add(u)
                    # DFS to find the connected components
                    while stack:
                        x = stack.pop()
                        for y in adj[x]:
                            if y in nbrs and y not in seen:
                                seen.add(y)
                                stack.append(y)
                τ[v] = comps
            return τ

        # Boundary 2: k > max_deg
        max_deg = max(len(adj[v]) for v in range(n))
        if k > max_deg:
            return τ       # All 0

        # ================== Normal Process ==================
        for v in range(n):
            nbrs = adj[v]
            if not nbrs:                     # Boundary 3: Isolated Vertex
                continue

            nbr_set = set(nbrs)

            # 1. Initialize local degree
            #Count the degree of each neighbor node u in the neighbor subgraph and store it in deg[u]
            deg = {u: 0 for u in nbrs}
            for u in nbrs:
                cnt = 0
                for w in adj[u]:
                    if w in nbr_set:
                        cnt += 1
                deg[u] = cnt

            # 2. k‑core peeling
            # Add neighbor nodes with degrees less than k to the stripping queue q;
            # Each time a node u is popped out from q, it is stripped from the graph (del deg[u])
            # and its neighbors' degrees are decremented by 1;
            # If a neighbor's degree becomes less than k, it is added to q;
            # This process continues until q is empty.
            q = deque([u for u, d in deg.items() if d < k])
            while q:
                u = q.popleft()
                if u not in deg:
                    continue
                for w in adj[u]:
                    if w in deg:
                        deg[w] -= 1
                        if deg[w] < k:
                            q.append(w)
                del deg[u]

            if not deg:                      # All deleted IF τ_k=0
                continue

            # 3. Count connected components 
            # Use DFS to calculate the number of connected components comps
            # Every time a new unvisited node u is encountered, it means it is a new connected component;
            # Then starting from this point, add all the nodes of this connected block to seen
            remaining = set(deg.keys())
            seen = set()
            comps = 0
            for u in remaining:
                if u in seen:
                    continue
                comps += 1
                stack = [u]
                seen.add(u)
                while stack:
                    x = stack.pop()
                    for y in adj[x]:
                        if y in remaining and y not in seen:
                            seen.add(y)
                            stack.append(y)
            τ[v] = comps                                     
            #For each node v, the number of connected components of its neighbor subgraph after k-core peeling is recorded in the result τ[v]  

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
