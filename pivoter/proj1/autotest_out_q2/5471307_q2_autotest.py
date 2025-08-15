#!/usr/bin/env python3
# Auto-generated for 5471307

STUDENT_ID = "5471307"
STUDENT_NAME = "Liangyu Li"

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
        # n = G.vertex_num
        # sd = [0] * n
        # return sd

        if hasattr(G, 'adj_list'):
            adj = G.adj_list
        elif hasattr(G, 'adj_list_out'):
            adj = G.adj_list_out
        else:
            raise AttributeError("Graph object must expose 'adj_list' or 'adj_list_out'.")
        n = G.vertex_num if hasattr(G, 'vertex_num') else len(adj)
        if n != len(adj):
            n = len(adj)
        tau = [0]*n
        if n == 0:
            return tau
        if k == 0:
            mark = [0]*n
            token = 0
            seen = [0]*n
            for v in range(n):
                token += 1
                for u in adj[v]:
                    mark[u] = token
                comp = 0
                for u in adj[v]:
                    if seen[u] != token:
                        comp += 1
                        stack = [u]
                        seen[u] = token
                        while stack:
                            x = stack.pop()
                            for w in adj[x]:
                                if mark[w] == token and seen[w] != token:
                                    seen[w] = token
                                    stack.append(w)
                tau[v] = comp
            return tau
        core = kCoreBaseStructuralDiversity._global_core_numbers(adj)
        mark = [0]*n
        alive = [False]*n
        int_deg = [0]*n
        removed = [False]*n
        token = 0
        for v in range(n):
            token += 1
            Sv = []
            for u in adj[v]:
                if core[u] < k:
                    continue
                mark[u] = token
                removed[u] = False
                Sv.append(u)
            if not Sv:
                tau[v] = 0
                continue
            for u in Sv:
                deg_in = 0
                for w in adj[u]:
                    if mark[w] == token:
                        deg_in += 1
                int_deg[u] = deg_in
                alive[u] = True
            dq = deque()
            for u in Sv:
                if int_deg[u] < k:
                    dq.append(u)
                    removed[u] = True
                    alive[u] = False
            while dq:
                u = dq.popleft()
                for w in adj[u]:
                    if mark[w] == token and not removed[w]:
                        int_deg[w] -= 1
                        if int_deg[w] == k-1:
                            removed[w] = True
                            alive[w] = False
                            dq.append(w)
            comp = 0
            for u in Sv:
                if alive[u]:
                    comp += 1
                    stack = [u]
                    alive[u] = False
                    while stack:
                        x = stack.pop()
                        for w in adj[x]:
                            if mark[w] == token and alive[w]:
                                alive[w] = False
                                stack.append(w)
            tau[v] = comp
        return tau

    @staticmethod
    def _global_core_numbers(adj):
        n = len(adj)
        deg = [len(adj[u]) for u in range(n)]
        if n == 0:
            return deg
        maxd = max(deg)
        bins = [0]*(maxd+1)
        for d in deg:
            bins[d] += 1
        start = 0
        for d in range(maxd+1):
            c = bins[d]
            bins[d] = start
            start += c
        vert = [0]*n
        pos = [0]*n
        for u in range(n):
            d = deg[u]
            idx = bins[d]
            pos[u] = idx
            vert[idx] = u
            bins[d] += 1
        for d in range(maxd, 0, -1):
            bins[d] = bins[d-1]
        bins[0] = 0
        for i in range(n):
            u = vert[i]
            for w in adj[u]:
                if deg[w] > deg[u]:
                    dw = deg[w]
                    pw = pos[w]
                    pw_start = bins[dw]
                    if w != vert[pw_start]:
                        w2 = vert[pw_start]
                        vert[pw_start], vert[pw] = vert[pw], vert[pw_start]
                        pos[w], pos[w2] = pw_start, pw
                    bins[dw] += 1
                    deg[w] -= 1
        return deg


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
