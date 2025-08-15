#!/usr/bin/env python3
# Auto-generated for 5519806

STUDENT_ID = "5519806"
STUDENT_NAME = "Jiaxin Xiong"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    _mark      = None  # int array length n, marks vertices that are neighbors of current v
    _kept_flag     = None # vertices that survive peeling
    _visited   = None     # BFS visited when counting components
    _deg_local = None     # degree inside neighbor induced subgraph
    _stamp     = 1        # vertex "version" across different iterations

    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Parameters: G, k
        Returns: List[int], tau_k(v) for all v
        """
        n   = G.vertex_num          # O(1)
        adj = G.adj_list            # O(1)

        #  each resize is O(n)
        K = kCoreBaseStructuralDiversity  # O(1)
        if K._mark is None or len(K._mark) < n:   # O(1) condition check
            # new array creations takes O(n) each
            K._mark      = [0] * n
            K._kept_flag     = [0] * n
            K._visited   = [0] * n
            K._deg_local = [0] * n


        tau = [0] * n   # O(n)

        # ---- Early eixt if k > maximum degree ---------------------------
        # max(...) loops over all vertices once, so O(n)
        if k > max((len(adj[v]) for v in range(n)), default=0):
            return tau  # O(1)

        # =================for each vertex v ======================
        # Loop runs n times; each iteration costs T_v.
        for v in range(n):                              # Condition check: n+1 times -> O(n)
            K._stamp += 1                               # O(1)
            stamp = K._stamp                            # O(1)

            nbrs = adj[v]                               # O(1)
            # If no neighbors, τ_k(v)=0 directly
            if not nbrs:                                # O(1)
                tau[v] = 0                              # O(1)
                continue                                # O(1)

            # -------- Mark neighbors of v -----------------------------
            # Loop d(v) times, for each iteration:
            #   condition check cost d(v)+1 -> O(d(v))
            for u in nbrs:                              # total O(d(v))
                K._mark[u] = stamp                      # O(1)

            # ---------Compute local degrees inside N(v) ----------------

            if k == 0:                                  # O(1)
                # Loop d(v) times, for each iteration:
                # condition check cost d(v)+1 -> O(d(v))
                for u in nbrs:                          # O(d(v))
                    K._kept_flag[u] = stamp                 # O(1)
                # Total O(d(v)).
            else:
                q = deque()                             # O(1)

                # 'For each neighbor u' loops d(v) times:
                #  iterate over adj[u] -> d(u)
                # Worst case cost: \sum_{u in N(v)} d(u)
                for u in nbrs:                          # loop check d(v)+1 -> O(d(v))
                    d = 0                               # O(1)
                    for w in adj[u]:                    # \sum d(u) iterations
                        # inner if: O(1)
                        if K._mark[w] == stamp:         # O(1)
                            d += 1                      # O(1)
                    K._deg_local[u] = d                 # O(1)
                    K._kept_flag[u] = stamp                 # O(1)
                    # if d<k push into queue: O(1)
                    if d < k:                           # O(1)
                        q.append(u)                     #  O(1)

                # ------------- Remove vertices with degree < k --------------
                # Each u that gets removed and then need to change the dgree for its neighbors.
                # A vertex can be enqueued multiple times but removed once.
                # Total edges processed <= |E_v|,
                # Total: O(\sum_{u \in N(v)} d(u))
                while q:                                # loop check t+1 times
                    u = q.popleft()                     # O(1)
                    if K._kept_flag[u] != stamp:        # O(1)
                        continue                        # O(1)
                    K._kept_flag[u] = 0                    # O(1)
                    for w in adj[u]:                    # d(u) iterations
                        if K._mark[w] == stamp and K._kept_flag[w] == stamp: #O(1)
                            K._deg_local[w] -= 1        # O(1)
                            if K._deg_local[w] == k - 1:
                                q.append(w)             # O(1)

            # ------------- Count connected components -------------
            # BFS over remaining kept vertices. Each vertex visited once.
            cnt = 0                                     # O(1)
            for u in nbrs:                              # O(d(v)) loop
                if K._kept_flag[u] == stamp and K._visited[u] != stamp:  # O(1)
                    cnt += 1                            # O(1)
                    K._visited[u] = stamp               # O(1)
                    bfs = deque([u])                    # O(1)
                    # BFS explores edges only between kept marked vertices.
                    # Total cost ≤ O(d(v) + |E_v|)
                    while bfs:                          # check visited vertices+1 -> O(|V_subg|)
                        x = bfs.popleft()               # O(1)
                        for w in adj[x]:                # d(x) iterations
                            if (K._mark[w] == stamp and
                                K._kept_flag[w] == stamp and
                                K._visited[w] != stamp):    # all O(1)
                                K._visited[w] = stamp       # O(1)
                                bfs.append(w)               # O(1)

            tau[v] = cnt                                # O(1)

        return tau   # O(1)

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
