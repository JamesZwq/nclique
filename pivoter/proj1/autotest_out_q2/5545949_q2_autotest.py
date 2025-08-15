#!/usr/bin/env python3
# Auto-generated for 5545949

STUDENT_ID = "5545949"
STUDENT_NAME = "Yuze Zou"

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
            induced_nodes = list(neighbors)
            induced_graph = kCoreBaseStructuralDiversity.induced_subgraph(G, induced_nodes)
            candidates = kCoreBaseStructuralDiversity.find_kcore(induced_graph, k)
            count = kCoreBaseStructuralDiversity.count_induced_graph(induced_graph, candidates)
            sd[v] = count

        return sd

    @staticmethod
    def induced_subgraph(G, nodes):

        node_map = {node:i for i, node in enumerate(nodes)}
        subgraph = [[] for _ in nodes]
        added = set()

        for u in nodes:
            for v in G.adj_list[u]:
                if v in node_map:
                    i, j = node_map[u], node_map[v]
                    if (i, j) not in added and (j, i) not in added:
                        subgraph[i].append(j)
                        subgraph[j].append(i)
                        added.add((i, j))
        return subgraph

    @staticmethod
    def find_kcore(adj_list, k):

        n = len(adj_list)
        degree = [len(neigh) for neigh in adj_list]

        max_deg = max(degree, default=0)
        bucket = [set() for _ in range(max_deg + 1)]
        for i in range(n):
            bucket[degree[i]].add(i)

        core = degree[:]
        upper = min(k, len(bucket))
        for d in range(upper):
            while bucket[d]:
                u = bucket[d].pop()
                for v in adj_list[u]:
                    if core[v] > d:
                        bucket[core[v]].discard(v)
                        core[v] -= 1
                        bucket[core[v]].add(v)

        return [core[i] >= k for i in range(n)]
    # def find_kcore(adj_list, k):

    #     n = len(adj_list)
    #     degrees = [len(adj_list[i]) for i in range(n)]
    #     visited = [False] * n
    #     queue = deque()

    #     for i in range(n):
    #         if degrees[i] < k:
    #             queue.append(i)
    #             visited[i] = True

    #     while queue:
    #         u = queue.popleft()
    #         for v in adj_list[u]:
    #             if not visited[v]:
    #                 degrees[v] -= 1
    #                 if degrees[v] < k:
    #                     queue.append(v)
    #                     visited[v] = True

    #     return [not visited[i] for i in range(n)]


    @staticmethod
    def count_induced_graph(adj_list, candidates):

        n = len(adj_list)
        visited = [False] * n
        count = 0

        for i in range(n):
            if not candidates[i] or visited[i]:
                continue
            count += 1
            queue = deque([i])
            visited[i] = True
            while queue:
                u = queue.popleft()
                for v in adj_list[u]:
                    if candidates[v] and not visited[v]:
                        visited[v] = True
                        queue.append(v)
        return count

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
