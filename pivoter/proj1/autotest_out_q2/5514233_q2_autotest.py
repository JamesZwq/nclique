#!/usr/bin/env python3
# Auto-generated for 5514233

STUDENT_ID = "5514233"
STUDENT_NAME = "(FUFU) Guangping Ren"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        adj = [list(dict.fromkeys(neigh)) for neigh in G.adj_list]  
        res = [0] * n

        for v in range(n):
            nbrs = list(dict.fromkeys(adj[v]))
            if len(nbrs) < k:
                continue

            sub_adj = kCoreBaseStructuralDiversity._subgraph(nbrs, adj)
            alive = kCoreBaseStructuralDiversity._peel(sub_adj, k)
            res[v] = kCoreBaseStructuralDiversity._ccount(sub_adj, alive)

        return res

    @staticmethod
    def _subgraph(nbrs, adj):
        """Constructs a neighbor-induced subgraph and returns an adjacency list under the local index"""
        m = len(nbrs)
        idx = {u: i for i, u in enumerate(nbrs)}
        nbrs_set = set(nbrs)
        sub_adj = [[] for _ in range(m)]

        for i, u in enumerate(nbrs):
            sub_adj[i] = list(dict.fromkeys(
                idx[w] for w in adj[u] if w in nbrs_set and w != u
            ))

        return sub_adj

    @staticmethod
    def _peel(sub_adj, k):
        """k-core peeled, returns whether it is alive"""
        m = len(sub_adj)
        deg = [len(nei) for nei in sub_adj]
        alive = [True] * m
        queue = deque(i for i, d in enumerate(deg) if d < k)

        while queue:
            u = queue.popleft()
            if not alive[u]:
                continue
            alive[u] = False
            for v in sub_adj[u]:
                if alive[v]:
                    deg[v] -= 1
                    if deg[v] < k:
                        queue.append(v)

        return alive

    @staticmethod
    def _ccount(sub_adj, alive):
        """Count the number of connected components of live nodes"""
        m = len(sub_adj)
        vis = [False] * m
        count = 0

        for i in range(m):
            if alive[i] and not vis[i]:
                count += 1
                stack = [i]
                vis[i] = True
                while stack:
                    u = stack.pop()
                    for v in sub_adj[u]:
                        if alive[v] and not vis[v]:
                            vis[v] = True
                            stack.append(v)

        return count

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
