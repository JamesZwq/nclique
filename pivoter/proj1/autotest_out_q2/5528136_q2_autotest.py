#!/usr/bin/env python3
# Auto-generated for 5528136

STUDENT_ID = "5528136"
STUDENT_NAME = "Pei Pei"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass


    # the neighbor subgraph G[N(v)] of a given node v
    @staticmethod
    def build_neighbor_subgraph(G, v):
        sub = set(G.adj_list[v])   # neighbors of v
        sub_adj = {u: [] for u in sub}
        deg = {}
        for u in sub:
            # Only keep neighbors within the neighbor set
            nbrs_in_sub = [w for w in G.adj_list[u] if w in sub]
            sub_adj[u] = nbrs_in_sub
            deg[u] = len(nbrs_in_sub)   # degree in the subgraph
        return sub_adj, deg


    # k-core peeling on the neighbor subgraph
    @staticmethod
    def peel_kcore(sub_adj, deg, k):
        # Initialize the queue with nodes of degree < k
        q_list = []
        for u in deg:
            d = deg[u]  # u's degree
            if d < k:
                q_list.append(u)  # <k, add to the list
        q = deque(q_list)
        removed = set()

        # Remove nodes
        while q:
            u = q.popleft()
            if u in removed:
                continue
            removed.add(u)
            # Decrease degrees of neighbors
            for w in sub_adj[u]:
                if w not in removed:
                    deg[w] -= 1
                    if deg[w] < k:
                        q.append(w)
        return set(deg) - removed  # remain


    # the number of connected components in the remaining k-core subgraph
    @staticmethod
    def count_components(remain, sub_adj):
        component, seen = 0, set()
        for u in remain:
            if u in seen:
                continue
            component += 1
            stack = [u]
            seen.add(u)
            # dfs
            while stack:
                x = stack.pop()
                for y in sub_adj[x]:
                    if y in remain and y not in seen:
                        seen.add(y)
                        stack.append(y)
        return component


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

        for v in range(n):
            # single node, neighbor graph is empty
            if not G.adj_list[v]:
                continue

            sub_adj, deg = kCoreBaseStructuralDiversity.build_neighbor_subgraph(G, v)
            remain = kCoreBaseStructuralDiversity.peel_kcore(sub_adj, deg, k)
            sd[v] = kCoreBaseStructuralDiversity.count_components(remain, sub_adj)

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
