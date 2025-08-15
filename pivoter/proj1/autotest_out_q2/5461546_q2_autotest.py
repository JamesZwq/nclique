#!/usr/bin/env python3
# Auto-generated for 5461546

STUDENT_ID = "5461546"
STUDENT_NAME = "Tung Shing Lai"

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

        n = G.vertex_num
        sd = [0] * n
        adj = G.adj_list

        for v in range(n):
            nbrs = adj[v]

            #map neighbour to local index
            idx = dict(zip(nbrs, range(len(nbrs))))

            #build adjacency among neighbours
            m = len(nbrs)
            nbr_adj = [[] for _ in range(m)]
            for i, u in enumerate(nbrs):
                for w in adj[u]:
                    if w in idx:
                        nbr_adj[i].append(idx[w])


            for i, nbr_list in enumerate(nbr_adj):
                orig = nbrs[i]
                linked = []
                for j in nbr_list:
                    linked.append(nbrs[j])
            #compute degree list in the induced subgraph
            deg = []
            for neighs in nbr_adj:
                deg.append(len(neighs))

            removed = kCoreBaseStructuralDiversity.peel_and_mark(deg, nbr_adj, k)
            #print(f"[DEBUG] removed mask: {removed}")
            comp_count = kCoreBaseStructuralDiversity.count_components(removed,nbr_adj)
            #assign structural diversity score to vertex v
            sd[v] = comp_count

        return sd


    @staticmethod
    def peel_and_mark(deg, nbr_adj, k):
        """
        Perform k-core peeling on the induced subgraph.
        vertex with degree < k are removed iteratively.
        """
        removed = [False] * len(deg)
        stack = []
        for i in range(len(deg)):
            if deg[i] < k:
                stack.append(i)

        while stack:
            u = stack.pop()
            if removed[u]:
                continue
            removed[u] = True
            for w in nbr_adj[u]:
                if not removed[w]:
                    deg[w] -= 1
                    if deg[w] < k:
                        stack.append(w)

        return removed

    @staticmethod
    def count_components(removed, nbr_adj):

        """
        Count the number of connected components in the remaining (non-removed) induced subgraph.
        Uses DFS.
        """

        visit = [False] * len(removed)
        comp = 0

        for i in range(len(removed)):
            if removed[i] or visit[i]:
                continue

            comp+=1
            stack = [i]
            visit[i] = True
            while stack:
                u = stack.pop()
                for we in nbr_adj[u]:
                    if not removed[we] and not visit[we]:
                        visit[we] = True
                        stack.append(we)
        return comp

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
