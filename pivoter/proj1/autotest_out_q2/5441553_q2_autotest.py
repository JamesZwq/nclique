#!/usr/bin/env python3
# Auto-generated for 5441553

STUDENT_ID = "5441553"
STUDENT_NAME = "Kayla Ridjab"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def find(parent, x):
        if parent[x] != x:
            parent[x] = kCoreBaseStructuralDiversity.find(parent, parent[x])
        return parent[x]

    @staticmethod
    def union(parent, size, x, y):
        x = kCoreBaseStructuralDiversity.find(parent, x)
        y = kCoreBaseStructuralDiversity.find(parent, y)
        if x == y:
            return
        if size[x] < size[y]:
            size[y] += size[x]
            parent[x] = y
        else:
            size[x] += size[y]
            parent[y] = x

    @staticmethod
    def ds_cc(parent, size, edges, vertices):
        for edge in edges:
            kCoreBaseStructuralDiversity.union(parent, size, edge[0], edge[1])

        connected_components = {}

        # create connected components sets/arrays
        for i in vertices:
          # point each vertex to itself
          connected_components[i] = [i]

        for i in vertices:
          if parent[i] != i:
            connected_components[kCoreBaseStructuralDiversity.find(parent, parent[i])].append(i)
            del connected_components[i]

        return connected_components.values()

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
        sd = [0] * n

        curr_v = 0

        for curr_v in range(0, G.vertex_num):

          edges = []
          for u in G.adj_list[curr_v]:
            for v in G.adj_list[u]:
              if v == curr_v or v not in G.adj_list[curr_v]:
                continue
              if sorted([u, v]) not in edges:
                edges.append(sorted([u, v]))

          # using disjoint set to get connected components
          # disjoint set code modified from given lecture code
          parent = [i for i in range(n)] # each vertex points to itself
          size = [1 for i in range(n)]

          ds_cc = kCoreBaseStructuralDiversity.ds_cc(parent, size, edges, G.adj_list[curr_v])

          # get τ_k(v)
          # filter out items in ds_cc whose length < k
          sd[curr_v] = len([i for i in ds_cc if len(i) > k])

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
