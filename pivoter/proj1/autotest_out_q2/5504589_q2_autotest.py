#!/usr/bin/env python3
# Auto-generated for 5504589

STUDENT_ID = "5504589"
STUDENT_NAME = "Jilei Zhu"

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
        res = [0] * n
        if k < 0:
            return res
        mk = [0] * n
        rm = [0] * n
        vis = [0] * n
        vm = 1
        vr = 1
        vv = 1
        deg = {}
        for v in range(n):
            Nv = adj[v]
            if len(Nv) <= k:
                res[v] = 0
                continue
            vm += 1
            for u in Nv:
                mk[u] = vm
            deg.clear()
            for u in Nv:
                d = 0
                for w in adj[u]:
                    if mk[w] == vm:
                        d += 1
                deg[u] = d
            vr += 1
            q = deque()
            for u in Nv:
                if deg[u] < k:
                    rm[u] = vr
                    q.append(u)
            while q:
                u = q.popleft()
                for w in adj[u]:
                    if mk[w] == vm and rm[w] != vr:
                        deg[w] -= 1
                        if deg[w] < k:
                            rm[w] = vr
                            q.append(w)
            vv += 1
            c = 0
            for u in Nv:
                if mk[u] == vm and rm[u] != vr and vis[u] != vv:
                    c += 1
                    dq = deque([u])
                    vis[u] = vv
                    while dq:
                        x = dq.popleft()
                        for w in adj[x]:
                            if mk[w] == vm and rm[w] != vr and vis[w] != vv:
                                vis[w] = vv
                                dq.append(w)
            res[v] = c
        return res


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
