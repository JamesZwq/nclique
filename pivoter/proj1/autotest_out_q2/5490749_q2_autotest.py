#!/usr/bin/env python3
# Auto-generated for 5490749

STUDENT_ID = "5490749"
STUDENT_NAME = "Shawn Wang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num             # Number of vertices in the graph
        adj = G.adj_list             # Adjacency list representation of the graph
        result = [0] * n             # Result array to store τ_k(v) for each node

        for v in range(n):
            neighbors = list(set(adj[v]))  # Neighbor set N(v), removing duplicates
            if not neighbors:
                result[v] = 0
                continue

            # Build induced subgraph G[N(v)]
            # Map global node IDs to local indices for the subgraph
            idx_map = {node: idx for idx, node in enumerate(neighbors)}
            rev_map = {idx: node for node, idx in idx_map.items()}
            size = len(neighbors)

            # Build adjacency list for the induced subgraph
            sub_adj = [[] for _ in range(size)]
            for i, u in enumerate(neighbors):
                for w in adj[u]:
                    if w in idx_map:
                        sub_adj[i].append(idx_map[w])

            # Perform k-core decomposition on the induced subgraph
            deg = [len(sub_adj[i]) for i in range(size)]  # Degree of each node in subgraph
            is_core = [True] * size                       # Mark nodes that are still in k-core
            queue = deque([i for i in range(size) if deg[i] < k])  # Nodes to remove initially
            for i in queue:
                is_core[i] = False

            # Iteratively remove nodes with degree < k
            while queue:
                u = queue.popleft()
                for w in sub_adj[u]:
                    if is_core[w]:
                        deg[w] -= 1
                        if deg[w] < k:
                            is_core[w] = False
                            queue.append(w)

            # Count connected components in the remaining k-core subgraph
            visited = [False] * size

            def bfs(start):
                q = deque([start])
                visited[start] = True
                while q:
                    u = q.popleft()
                    for w in sub_adj[u]:
                        if is_core[w] and not visited[w]:
                            visited[w] = True
                            q.append(w)

            count = 0
            for i in range(size):
                if is_core[i] and not visited[i]:
                    bfs(i)
                    count += 1

            result[v] = count  # τ_k(v) = number of connected components in k-core of G[N(v)]

        return result

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
