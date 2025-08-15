#!/usr/bin/env python3
# Auto-generated for 5529212

STUDENT_ID = "5529212"
STUDENT_NAME = "Chen Lei"

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
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            neighbors = set(neighbors)
            if not neighbors:
                sd[v] = 0
                continue

            # create Neighbour-induced subgraph of v
            neighbours_G = {u:[] for u in neighbors}
            for u in neighbors:
                for neighbour_u in G.adj_list[u]:
                    if neighbour_u in neighbors:
                        neighbours_G[u].append(neighbour_u)

            sd[v] = kCoreBaseStructuralDiversity.kcore(neighbours_G, k)
        return sd


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def kcore(G, k):
        if not G:
            return 0
        
        deg = {u:len(neighbors) for u, neighbors in G.items()}
        deleted_nodes = set()
        updated_nodes = deque()
        
        for v in G:
            if deg[v] < k:
                deleted_nodes.add(v)
                updated_nodes.append(v)


        while updated_nodes:
            u = updated_nodes.popleft()
            for v in G[u]:
                if v not in deleted_nodes:
                    deg[v] -= 1
                    if deg[v] < k:
                        deleted_nodes.add(v)
                        updated_nodes.append(v)

        remains = [v for v in G if v not in deleted_nodes]
        if not remains:
            return 0
        
        # count k_core by bfs
        visited = set()
        kcore_cnt = 0

        for node in remains:
            if node not in visited:
                kcore_cnt += 1
                q = deque()
                q.append(node)

                while q:
                    u = q.popleft()
                    visited.add(u)

                    for v in G[u]:
                        if v not in deleted_nodes and v not in visited:
                            q.append(v)
        return kcore_cnt

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
