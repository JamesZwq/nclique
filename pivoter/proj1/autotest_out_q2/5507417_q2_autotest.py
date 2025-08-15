#!/usr/bin/env python3
# Auto-generated for 5507417

STUDENT_ID = "5507417"
STUDENT_NAME = "Chenqi Zhang"

# ======= 学生代码 =======
from collections import deque
from typing import List, Set


class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k: int) -> List[int]:
        n = G.vertex_num
        result = [0] * n
        for v in range(n):
            # Get all neighbors of node v
            nbrs = G.adj_list[v]
            subgraph_nodes = set(nbrs)

            # Build subgraph induced by neighbors of v
            subgraph = kCoreBaseStructuralDiversity.build_subgraph(G, subgraph_nodes)

            # Get all nodes that are in k-core
            k_core_nodes = kCoreBaseStructuralDiversity.find_k_core(subgraph, k)

            # Find connected components among k-core nodes
            components = kCoreBaseStructuralDiversity.find_connected_components(k_core_nodes, subgraph)

            # The number of k-core components is the diversity of node v
            result[v] = len(components)

        return result

    @staticmethod
    def build_subgraph(G, nodes: Set[int]) -> List[List[int]]:
        subgraph = [[] for _ in range(G.vertex_num)]
        for u in nodes:
            for v in G.adj_list[u]:
                if v in nodes:
                    subgraph[u].append(v)
        return subgraph

    @staticmethod
    def find_k_core(adj_list: List[List[int]], k: int) -> Set[int]:
        degree = [len(neigh) for neigh in adj_list]
        q = deque([i for i in range(len(adj_list)) if degree[i] < k])
        removed = set()

        while q:
            u = q.popleft()
            if u in removed:
                continue
            removed.add(u)
            for v in adj_list[u]:
                degree[v] -= 1
                if degree[v] == k - 1:
                    q.append(v)

        result = set()
        for i in range(len(adj_list)):
            if i not in removed:
                if len(adj_list[i]) >= k:
                    result.add(i)
        return result



    @staticmethod
    def find_connected_components(nodes: Set[int], adj_list: List[List[int]]) -> List[Set[int]]:
        visited = set()
        components = []

        def bfs(start):
            q = deque([start])
            comp = {start}
            visited.add(start)
            while q:
                u = q.popleft()
                for v in adj_list[u]:
                    if v in nodes and v not in visited:
                        visited.add(v)
                        comp.add(v)
                        q.append(v)
            return comp

        for v in nodes:
            if v not in visited:
                components.append(bfs(v))
        return components

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
