#!/usr/bin/env python3
# Auto-generated for 5418256

STUDENT_ID = "5418256"
STUDENT_NAME = "Sai Nair"

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
        sd = []

        # Calculate degree for each vertex in neighbour induced subgraph for vertices
        result = [[0 for _ in range(n)] for _ in range(n)]
        neighbours = [set(G.adj_list[v]) for v in range(n)]
        for v in range(n):
            for u in neighbours[v]:
                if u < v: continue
                degree_in_induced = len(neighbours[v] & neighbours[u])
                result[v][u] = degree_in_induced
                result[u][v] = degree_in_induced

        # Calculate k-core structural diversity for each vertex
        for i in range(n):
            # Initialise neighbours as in the k-core
            in_core = [False] * n
            for j in G.adj_list[i]:
                in_core[j] = True

            # Calculate k-core of neighbour induced subgraph
            deg = result[i]
            stack = [u for u in range(n) if deg[u] < k]
            while stack:
                u = stack.pop()
                if not in_core[u]: continue
                in_core[u] = False
                for v in G.adj_list[u]:
                    if in_core[v]:
                        deg[v] -= 1
                        if deg[v] == k - 1:
                            stack.append(v)

            # Count connected components in k-core
            # in_core now acts as a negation of the visited array - in_core[u] = True means u hasn't been visited yet
            count = 0
            for u in range(n):
                if not in_core[u]: continue
                count += 1
                in_core[u] = False
                stack = [u]
                while stack:
                    cur = stack.pop()
                    for v in G.adj_list[cur]:
                        if not in_core[v]: continue
                        in_core[v] = False
                        stack.append(v)
            sd.append(count)

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
