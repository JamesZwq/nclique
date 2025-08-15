#!/usr/bin/env python3
# Auto-generated for 5560661

STUDENT_ID = "5560661"
STUDENT_NAME = "Cesare Li"

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
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        # TODO
        τ = [0] * G.vertex_num
        for v in range(G.vertex_num):
            # Step 1: Get neighbors of v
            neighbors = set(G.adj_list[v])
            if not neighbors:
                τ[v] = 0
                continue

            # Step 2: Initialize S as N(v) and compute degrees within S
            S = neighbors.copy()
            deg_S = {u: 0 for u in S}
            for u in S:
                for w in G.adj_list[u]:
                    if w in S:
                        deg_S[u] += 1

            # Step 3: Compute k-core by removing vertices with degree < k
            queue = deque()
            in_queue = set()
            for u in S:
                if deg_S[u] < k:
                    queue.append(u)
                    in_queue.add(u)

            while queue:
                u = queue.popleft()
                in_queue.remove(u)
                if u not in S:
                    continue
                S.remove(u)
                for w in G.adj_list[u]:
                    if w in S:
                        deg_S[w] -= 1
                        if deg_S[w] < k and w not in in_queue:
                            queue.append(w)
                            in_queue.add(w)

            # Step 4: Count connected components in the k-core
            if not S:
                τ[v] = 0
            else:
                τ[v] = kCoreBaseStructuralDiversity.count_connected_components(G, list(S))

        return τ

    @staticmethod
    def count_connected_components(G, vertices):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        vertices : List[int]  # Subset of vertices
        Returns
        -------
        int  # Number of connected components
        """
        parent = {v: v for v in vertices}
        rank = {v: 0 for v in vertices}

        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x, y):
            px, py = find(x), find(y)
            if px != py:
                if rank[px] > rank[py]:
                    parent[py] = px
                elif rank[px] < rank[py]:
                    parent[px] = py
                else:
                    parent[py] = px
                    rank[px] += 1

        for v in vertices:
            for u in G.adj_list[v]:
                if u in vertices:
                    union(v, u)

        components = set(find(v) for v in vertices)
        return len(components)

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
