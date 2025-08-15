#!/usr/bin/env python3
# Auto-generated for 5514565

STUDENT_ID = "5514565"
STUDENT_NAME = "Tiying Song"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from typing import List, Dict
################################################################################

class kCoreBaseStructuralDiversity(object):
    # ------------------------------------------------------------------ #
    # Helper functions
    # ------------------------------------------------------------------ #
    @staticmethod
    def _build_local_id(nbrs: List[int]) -> Dict[int, int]:
        return {gid: idx for idx, gid in enumerate(nbrs)}

    @staticmethod
    def _k_core_components(adj_loc: List[List[int]], k: int) -> int:
        d = len(adj_loc)
        if d == 0:
            return 0

        # ---- k‑core peeling (O(m_sub)) ----
        deg = [len(nei) for nei in adj_loc]
        in_core = [True] * d
        stack = [i for i, dv in enumerate(deg) if dv < k]
        while stack:
            u = stack.pop()
            if not in_core[u]:
                continue
            in_core[u] = False
            for v in adj_loc[u]:
                if in_core[v]:
                    deg[v] -= 1
                    if deg[v] < k:
                        stack.append(v)

        # ---- count connected components on the remaining vertices ----
        seen = [False] * d
        cc = 0
        for i in range(d):
            if in_core[i] and not seen[i]:
                cc += 1
                q = deque([i])
                seen[i] = True
                while q:
                    u = q.popleft()
                    for v in adj_loc[u]:
                        if in_core[v] and not seen[v]:
                            seen[v] = True
                            q.append(v)
        return cc

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #
    @staticmethod
    def process(G, k: int) -> List[int]:
        n = G.vertex_num
        sd = [0] * n
        adj_global = G.adj_list

        for v in range(n):
            nbrs = adj_global[v]
            d = len(nbrs)

            if d < k:
                sd[v] = 0
                continue

            # ---- Build neighbour
            local_id = kCoreBaseStructuralDiversity._build_local_id(nbrs)
            adj_loc = [[] for _ in range(d)]
            for g_u in nbrs:
                for g_v in adj_global[g_u]:
                    if g_v in local_id:
                        u = local_id[g_u]
                        w = local_id[g_v]
                        if w not in adj_loc[u]:
                            adj_loc[u].append(w)
                        if u not in adj_loc[w]:
                            adj_loc[w].append(u)

            # ---- Compute τ_k(v) on the local graph ----
            sd[v] = kCoreBaseStructuralDiversity._k_core_components(adj_loc, k)

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
