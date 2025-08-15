#!/usr/bin/env python3
# Auto-generated for 5504665

STUDENT_ID = "5504665"
STUDENT_NAME = "Zilong Xia"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    @staticmethod
    def process(G, k):
        """
        Compute τₖ(v) for every vertex v in G, where τₖ(v) is the number of connected components
        in the k-core of v's ego‑network.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            Graph with attributes:
                vertex_num : int
                adj_list   : List[List[int]]
        k : int
            Core threshold.

        Returns
        -------
        List[int]
            τₖ(v) for each v in range(G.vertex_num).
        """
        n = G.vertex_num
        adj = G.adj_list
        sd = [0] * n

        # Process each vertex independently (embarrassingly parallelizable)
        for v in range(n):
            neighbors = adj[v]
            # If v's degree < k, its ego‑network k-core is empty → τₖ(v)=0
            if len(neighbors) < k:
                continue

            # Build index mapping for neighbors of v
            U = neighbors
            idx = {u: i for i, u in enumerate(U)}
            m = len(U)

            # Build induced subgraph on U
            degree = [0] * m
            neighbors_idx = [[] for _ in range(m)]
            for i, u in enumerate(U):
                for w in adj[u]:
                    j = idx.get(w)
                    if j is not None:
                        neighbors_idx[i].append(j)
                degree[i] = len(neighbors_idx[i])

            # Peeling to local k-core
            removed = [False] * m
            queue = deque([i for i, d in enumerate(degree) if d < k])
            while queue:
                i = queue.popleft()
                if removed[i]:
                    continue
                removed[i] = True
                for j in neighbors_idx[i]:
                    if not removed[j]:
                        degree[j] -= 1
                        if degree[j] == k - 1:
                            queue.append(j)

            # Count connected components among remaining nodes
            component_id = [-1] * m
            comp_index = 0
            for i in range(m):
                if not removed[i] and component_id[i] == -1:
                    stack = [i]
                    component_id[i] = comp_index
                    while stack:
                        x = stack.pop()
                        for y in neighbors_idx[x]:
                            if not removed[y] and component_id[y] == -1:
                                component_id[y] = comp_index
                                stack.append(y)
                    comp_index += 1

            sd[v] = comp_index

        return sd

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
