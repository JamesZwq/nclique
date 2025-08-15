#!/usr/bin/env python3
# Auto-generated for 5520900

STUDENT_ID = "5520900"
STUDENT_NAME = "Yucan Liu"

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
        Compute τₖ(v) for every vertex v in an undirected, unweighted graph G.
        G.vertex_num gives |V|, and G.adj_list is a list of neighbor lists.
        """
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neigh = G.adj_list[v]
            dv = len(neigh)
            # If v has no neighbors, its k‑core is empty → diversity = 0
            if dv == 0:
                sd[v] = 0
                continue

            # Build an index for the induced subgraph on N(v)
            idx = {u: i for i, u in enumerate(neigh)}
            induced_adj = [[] for _ in range(dv)]
            # Populate the induced‐subgraph adjacency
            for u in neigh:
                ui = idx[u]
                for w in G.adj_list[u]:
                    if w in idx:               # only keep neighbors in N(v)
                        induced_adj[ui].append(idx[w])

            # Peeling process to compute the k‑core of this induced subgraph
            deg = [len(lst) for lst in induced_adj]
            removed = [False] * dv
            queue = deque(i for i, d in enumerate(deg) if d < k)
            while queue:
                i = queue.popleft()
                if removed[i]:
                    continue
                removed[i] = True
                for j in induced_adj[i]:
                    if not removed[j]:
                        deg[j] -= 1
                        if deg[j] == k - 1:
                            queue.append(j)

            # Gather the surviving core nodes
            core_nodes = [i for i in range(dv) if not removed[i]]
            if not core_nodes:
                sd[v] = 0
                continue

            # Count connected components in the k‑core
            visited = [False] * dv
            comps = 0
            for i in core_nodes:
                if not visited[i]:
                    comps += 1
                    dq = deque([i])
                    visited[i] = True
                    while dq:
                        x = dq.popleft()
                        for j in induced_adj[x]:
                            if not removed[j] and not visited[j]:
                                visited[j] = True
                                dq.append(j)

            sd[v] = comps

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
