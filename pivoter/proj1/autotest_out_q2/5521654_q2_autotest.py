#!/usr/bin/env python3
# Auto-generated for 5521654

STUDENT_ID = "5521654"
STUDENT_NAME = "Jiayou Wang"

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
        n = G.vertex_num
        tau = [0] * n  # τ_k(v)

        for v in range(n):
            nbrs = G.adj_list[v]          # N(v)

            # Quick prune: impossible to have a k-core
            if len(nbrs) < k:
                tau[v] = 0
                continue

            # ---- Build neighbour-induced subgraph on N(v) (exclude v itself) ----
            idx = {u: i for i, u in enumerate(nbrs)}  # original id -> [0..m-1]
            m = len(nbrs)
            sub_adj = [[] for _ in range(m)]
            nbrs_set = set(nbrs)

            for u in nbrs:
                ui = idx[u]
                for w in G.adj_list[u]:
                    if w in nbrs_set:
                        sub_adj[ui].append(idx[w])

            # Count the number of connected components in the k-core
            tau[v] = kCoreBaseStructuralDiversity._count_k_core(sub_adj, k)

        return tau

    @staticmethod
    def _count_k_core(adj, k):
        """
        Return the number of connected components inside the k-core of `adj`.
        `adj` is a 0..m-1 indexed adjacency list.
        """
        m = len(adj)
        if m == 0:
            return 0

        # 1) k-core peeling
        deg = [len(adj[u]) for u in range(m)]
        pruned = [False] * m
        q = deque([u for u in range(m) if deg[u] < k])
        for u in q:
            pruned[u] = True

        while q:
            u = q.popleft()
            for v in adj[u]:
                if not pruned[v]:
                    deg[v] -= 1
                    if deg[v] < k:
                        pruned[v] = True
                        q.append(v)

        core_nodes = [u for u in range(m) if not pruned[u]]
        if not core_nodes:
            return 0

        # 2) Count components within the remaining k-core
        seen = [False] * m
        comp_cnt = 0
        for s in core_nodes:
            if not seen[s]:
                comp_cnt += 1
                dq = deque([s])
                seen[s] = True
                while dq:
                    u = dq.popleft()
                    for v in adj[u]:
                        if not pruned[v] and not seen[v]:
                            seen[v] = True
                            dq.append(v)

        return comp_cnt



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
