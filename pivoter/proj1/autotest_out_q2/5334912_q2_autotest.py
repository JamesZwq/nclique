#!/usr/bin/env python3
# Auto-generated for 5334912

STUDENT_ID = "5334912"
STUDENT_NAME = "Zhihan Sun"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        # Iterate through each vertex v, get its neighbors from the adjacency list.
        for v in range(n):
            nbrs = G.adj_list[v]
            if not nbrs:
                sd[v] = 0
                continue
            else:
                # build the neighbor-induced subgraph and compute k-core components.
                nbrs = set(nbrs)
                sub_G = kCoreBaseStructuralDiversity._build_subgraph(G, nbrs)
                sd[v] = kCoreBaseStructuralDiversity.kcore_process(sub_G, k)

        return sd

    @staticmethod
    def kcore_process(sub_G, k):
        if not sub_G:
            return 0

        # Initialize degree of each node in the subgraph.
        degrees = {}
        for u in sub_G:
            length = len(sub_G[u])
            degrees[u] = length

        de_list = []

        vertices_all = deque()
        for u in sub_G:
            if degrees[u] < k:
                de_list.append(u)
                vertices_all.append(u)

        # remove nodes with degree < k and update neighbors' degrees.
        while vertices_all:
            u = vertices_all.popleft()
            for vertex in sub_G[u]:
                if vertex not in de_list:
                    degrees[vertex] -= 1
                    if k > degrees[vertex]:
                        de_list.append(vertex)
                        vertices_all.append(vertex)

        return kCoreBaseStructuralDiversity._process_kcore_components(sub_G, de_list)



    @staticmethod
    def _build_subgraph(G, nbrs):

        # Build neighbor-induced subgraph using given neighbor set nbrs.
        sub_G = {}
        for u in nbrs:
            sub_G[u] = []
        for u in nbrs:
            for u_nbr in G.adj_list[u]:
                if u_nbr in nbrs:
                    sub_G[u].append(u_nbr)
        return sub_G

    @staticmethod
    def _process_kcore_components(sub_G, de_list):

        # Gather remaining vertices that are still in the subgraph.
        remained_vertices = []
        for u in sub_G:
            if u not in de_list:
                remained_vertices.append(u)

        if not remained_vertices:
            return 0

        result = 0
        visited = []

        for node in remained_vertices:
            if node not in visited:
                result += 1
                kCoreBaseStructuralDiversity._dfs(node, sub_G, de_list, visited)

        return result

    @staticmethod
    def _dfs(start, sub_G, de_list, visited):
        
        # Use DFS from each unvisited node to count connected components.
        stack = []
        stack.append(start)
        while stack:
            u = stack.pop()
            if u in visited:
                continue
            visited.append(u)
            for vertex in sub_G[u]:
                if vertex not in de_list:
                    if vertex not in visited:
                        stack.append(vertex)

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
