#!/usr/bin/env python3
# Auto-generated for 5527851

STUDENT_ID = "5527851"
STUDENT_NAME = "Haoqian He"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    @staticmethod
    def process(G, k):
        """
        Global k-core prefiltering + local k-core peeling + connected component counting
        Inputs:
          - G: UndirectedUnweightedGraph, provides G.vertex_num and G.adj_list
          - k: the order of the k-core to compute
        Output:
          - tau: a list of length n, where tau[v] = number of connected components in the k-core
                 of the subgraph induced by v's neighbors
        """
        n = G.vertex_num
        adj = G.adj_list  # List[List[int]]

        # 1) Global k-core decomposition (bucket peeling), complexity: O(m + n)
        deg = [len(adj[u]) for u in range(n)]
        max_deg = max(deg) if n > 0 else 0
        buckets = [deque() for _ in range(max_deg + 1)]
        for u in range(n):
            buckets[deg[u]].append(u)

        core_number = [0] * n
        curr_deg = deg[:]  # dynamic degree array
        for d in range(max_deg + 1):
            while buckets[d]:
                u = buckets[d].popleft()
                core_number[u] = d
                for v in adj[u]:
                    if curr_deg[v] > d:
                        old = curr_deg[v]
                        curr_deg[v] -= 1
                        buckets[old].append(v)

        # 2) For each vertex v: filter neighbors with global core >= k, then perform local k-core peeling
        tau = [0] * n
        visited = [False] * n

        for v in range(n):
            # 2.1 Retain neighbors whose global core_number >= k
            Nk = [u for u in adj[v] if core_number[u] >= k]
            if not Nk:
                tau[v] = 0
                continue

            S = set(Nk)

            # 2.2 Local k-core peeling:
            #     First, compute local degree within subgraph
            deg_local = {}
            for u in Nk:
                cnt = 0
                for w in adj[u]:
                    if w in S:
                        cnt += 1
                deg_local[u] = cnt

            #     Remove vertices with degree < k in batches
            q = deque(u for u in Nk if deg_local[u] < k)
            while q:
                u = q.popleft()
                if u not in S:
                    continue
                S.remove(u)
                for w in adj[u]:
                    if w in S:
                        deg_local[w] -= 1
                        if deg_local[w] == k - 1:
                            q.append(w)

            if not S:
                tau[v] = 0
                continue

            # 2.3 On remaining nodes in S, perform BFS to count connected components
            comp_cnt = 0
            for u in S:
                if not visited[u]:
                    comp_cnt += 1
                    dq = deque([u])
                    visited[u] = True
                    while dq:
                        x = dq.popleft()
                        for w in adj[x]:
                            if w in S and not visited[w]:
                                visited[w] = True
                                dq.append(w)
            # Reset visited status
            for u in S:
                visited[u] = False

            tau[v] = comp_cnt

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
