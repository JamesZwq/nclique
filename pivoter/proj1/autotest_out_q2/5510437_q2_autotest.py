#!/usr/bin/env python3
# Auto-generated for 5510437

STUDENT_ID = "5510437"
STUDENT_NAME = "Yuhao Zhu"

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
        n = G.vertex_num
        if n == 0:
            return []
        # Create an adjacency list from the graph
        adj_list = [[] for _ in range(n)]
        # Reconstruct edges from the graph's representation
        for u in range(n):
            for v in G.adj_list[u]:
                if u < v:  # Avoid duplicate edges
                    adj_list[u].append(v)
                    adj_list[v].append(u)

        # compute structural diversity using the adjacency list
        sd = [0] * n
        mark = [False] * n
        deg_in_H = [0] * n

        for v in range(n):
            S = adj_list[v]
            d_v = len(S)
            if d_v == 0:
                sd[v] = 0
                continue
            if k > 0 and k > d_v - 1:
                sd[v] = 0
                continue

            # Mark all neighbors of v
            for u in S:
                mark[u] = True

            # Calculate degrees in the induced subgraph H
            for u in S:
                cnt = 0
                for w in adj_list[u]:
                    if mark[w]:
                        cnt += 1
                deg_in_H[u] = cnt

            # k-core decomposition using BFS
            remaining = set(S)
            q = deque()
            for u in S:
                if deg_in_H[u] < k:
                    q.append(u)

            while q:
                u = q.popleft()
                if u not in remaining:
                    continue
                remaining.remove(u)
                for w in adj_list[u]:
                    if mark[w] and w in remaining:
                        deg_in_H[w] -= 1
                        if deg_in_H[w] < k:
                            q.append(w)

            # Count connected components in the k-core
            visited = set()
            comp_count = 0
            for u in remaining:
                if u not in visited:
                    comp_count += 1
                    stack = [u]
                    visited.add(u)
                    while stack:
                        cur = stack.pop()
                        for w in adj_list[cur]:
                            if w in remaining and w not in visited:
                                visited.add(w)
                                stack.append(w)
            sd[v] = comp_count

            # Reset marks for next vertex
            for u in S:
                mark[u] = False

        return sd


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
