#!/usr/bin/env python3
# Auto-generated for 5534620

STUDENT_ID = "5534620"
STUDENT_NAME = "Haoqing Gu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from typing import List, Set, Dict

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
        sd = [0] * n
        if k < 0 or n == 0:
            return sd


        core_num = kCoreBaseStructuralDiversity._compute_core_numbers(G)
        if max(core_num) < k + 1:
            return sd

        k1_vertices = {v for v, c in enumerate(core_num) if c >= k + 1}
        comps = kCoreBaseStructuralDiversity._find_components(G, k1_vertices)

        for comp in comps:
            cand = set(comp)
            for v in comp:
                cand.update(G.adj_list[v])
            sd_local = kCoreBaseStructuralDiversity._process_component_batch(G, cand, comp, k)
            for v, val in sd_local.items():
                sd[v] = max(sd[v], val)
        return sd

    @staticmethod
    def _compute_core_numbers(G):
        n = G.vertex_num
        degree = [len(G.adj_list[v]) for v in range(n)]
        max_deg = max(degree) if degree else 0
        cnt = [0] * (max_deg + 1)
        for d in degree:
            cnt[d] += 1
        start = [0] * (max_deg + 1)
        acc = 0
        for d in range(max_deg + 1):
            start[d] = acc
            acc += cnt[d]
        vo = [0] * n
        pos = start.copy()
        for v, d in enumerate(degree):
            vo[pos[d]] = v
            pos[d] += 1
        pos = start.copy()
        for i in range(n):
            v = vo[i]
            for u in G.adj_list[v]:
                if degree[u] > degree[v]:
                    du = degree[u]
                    pu = pos[du]
                    w = vo[pu]
                    if u != w:
                        vo[pu] = u
                        vo[i] = w
                        pos[du] += 1
                    degree[u] -= 1
        return degree

    @staticmethod
    def _process_component_batch(G, cand, comp_core, k):
        sd_local = {}
        cache = {}
        for v in cand:
            nbr = set(G.adj_list[v])
            cache[v] = {u for u in nbr if len(G.adj_list[u]) >= k}
        for v in comp_core:
            nbr_set = cache[v]
            if len(nbr_set) < k:
                sd_local[v] = 0
            else:
                sd_local[v] = kCoreBaseStructuralDiversity._compute_k_core_count(G, nbr_set, cache, k)
        return sd_local

    @staticmethod
    def _compute_k_core_count(G, verts, cache, k):
        if len(verts) < k:
            return 0
        deg = {v: len(cache[v] & verts) for v in verts}
        q = deque([v for v, d in deg.items() if d < k])
        rem = set(q)
        while q:
            v = q.popleft()
            for u in cache[v] & verts:
                if u not in rem:
                    deg[u] -= 1
                    if deg[u] < k:
                        q.append(u)
                        rem.add(u)
        remain = verts - rem
        if not remain:
            return 0
        return kCoreBaseStructuralDiversity._count_components_set(G, remain)

    @staticmethod
    def _count_components_set(G, verts):
        seen = set()
        cnt = 0
        for v in verts:
            if v in seen:
                continue
            stack = [v]
            seen.add(v)
            while stack:
                u = stack.pop()
                for w in G.adj_list[u]:
                    if w in verts and w not in seen:
                        seen.add(w)
                        stack.append(w)
            cnt += 1
        return cnt

    @staticmethod
    def _find_components(G, verts):
        if not verts:
            return []
        vis = set()
        comps = []
        for v in verts:
            if v in vis:
                continue
            comp = set()
            dq = deque([v])
            vis.add(v)
            while dq:
                u = dq.popleft()
                comp.add(u)
                for w in G.adj_list[u]:
                    if w in verts and w not in vis:
                        vis.add(w)
                        dq.append(w)
            comps.append(comp)
        return comps

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
