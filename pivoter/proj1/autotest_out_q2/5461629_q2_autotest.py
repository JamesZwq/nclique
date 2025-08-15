#!/usr/bin/env python3
# Auto-generated for 5461629

STUDENT_ID = "5461629"
STUDENT_NAME = "Scarlett Luo"

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
        adj = G.adj_list
        if n == 0:
            return []

        core = kCoreBaseStructuralDiversity._core_numbers(adj)

        sd = [0] * n
        mark = [0] * n
        deg_local = [0] * n
        token = 1

        for v in range(n):
            cand = []
            for u in adj[v]:
                if core[u] >= k:
                    mark[u] = token
                    deg_local[u] = 0
                    cand.append(u)

            m = len(cand)
            if m < k:
                sd[v] = 0
                token += 1
                continue

            for u in cand:
                for w in adj[u]:
                    if mark[w] == token:
                        deg_local[u] += 1

            q = deque()
            for u in cand:
                if deg_local[u] < k:
                    q.append(u)

            dead_token = token + 1
            while q:
                u = q.popleft()
                if mark[u] != token:
                    continue
                mark[u] = dead_token
                for w in adj[u]:
                    if mark[w] == token:
                        deg_local[w] -= 1
                        if deg_local[w] == k - 1:
                            q.append(w)

            comp_token = token + 2
            components = 0
            for u in cand:
                if mark[u] == token:
                    components += 1
                    stack = [u]
                    mark[u] = comp_token
                    while stack:
                        x = stack.pop()
                        for w in adj[x]:
                            if mark[w] == token:
                                mark[w] = comp_token
                                stack.append(w)

            sd[v] = components
            token += 3

            if token > 1_000_000_000:
                mark = [0] * n
                token = 1

        return sd

    @staticmethod
    def _core_numbers(adj):
        n = len(adj)
        deg = [len(nei) for nei in adj]
        if n == 0:
            return []

        max_deg = max(deg)
        bins = [deque() for _ in range(max_deg + 1)]
        for v, d in enumerate(deg):
            bins[d].append(v)

        core = [-1] * n
        for d in range(max_deg + 1):
            while bins[d]:
                v = bins[d].pop()
                if core[v] != -1:
                    continue
                core[v] = d
                for u in adj[v]:
                    if core[u] == -1:
                        du = deg[u]
                        if du > 0:
                            deg[u] = du - 1
                            bins[du - 1].append(u)
        return core



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
