#!/usr/bin/env python3
# Auto-generated for 5349040

STUDENT_ID = "5349040"
STUDENT_NAME = "Ziqi Zeng"

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
        Compute τ_k(v) for every vertex v:
        For each v, build the neighbor-induced subgraph, peel nodes with degree < k,
        then count the number of remaining connected components.
        """
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            # ----- Build neighbor-induced subgraph (nodes & edges) -----
            neigh = G.adj_list[v]
            if len(neigh) == 0:
                sd[v] = 0
                continue

            nbr_set = set()
            for u in neigh:
                nbr_set.add(u)

            nbr_adj = {}
            for u in neigh:
                nbr_adj[u] = []

            for u in neigh:
                adj_u = G.adj_list[u]
                for w in adj_u:
                    in_set = (w in nbr_set)
                    is_center = (w == v)
                    if in_set and (not is_center):
                        nbr_adj[u].append(w)

            # ----- Peel nodes with degree < k (k-core pruning) -----
            if k > 0:
                q = deque()
                for u in nbr_adj:
                    deg_u = len(nbr_adj[u])
                    if deg_u < k:
                        q.append(u)

                removed = set()
                while len(q) > 0:
                    u = q.popleft()
                    if u in removed:
                        continue
                    removed.add(u)

                    if u in nbr_adj:
                        neighbors_of_u = nbr_adj[u]
                    else:
                        neighbors_of_u = []

                    for w in neighbors_of_u:
                        if w in removed:
                            continue
                        if w in nbr_adj:
                            list_w = nbr_adj[w]
                            if u in list_w:
                                list_w.remove(u)
                                deg_w = len(list_w)
                                if deg_w < k:
                                    q.append(w)

                    if u in nbr_adj:
                        del nbr_adj[u]

            # ----- Count connected components in the remaining subgraph -----
            comp_cnt = 0
            visited = set()

            keys_snapshot = list(nbr_adj.keys())
            for u in keys_snapshot:
                if u in visited:
                    continue

                comp_cnt += 1
                dq = deque()
                dq.append(u)
                visited.add(u)

                while len(dq) > 0:
                    x = dq.popleft()
                    for w in nbr_adj[x]:
                        exists_w = (w in nbr_adj)
                        not_visited = (w not in visited)
                        if exists_w and not_visited:
                            visited.add(w)
                            dq.append(w)

            sd[v] = comp_cnt

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
