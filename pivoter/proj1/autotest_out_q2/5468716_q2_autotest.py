#!/usr/bin/env python3
# Auto-generated for 5468716

STUDENT_ID = "5468716"
STUDENT_NAME = "Yiyang Qian"

# ======= 学生代码 =======
from collections import deque, defaultdict
from typing import List

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k: int) -> List[int]:
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The undirected, unweighted input graph.
        k : int
            The minimum core number to consider.

        Returns
        -------
        List[int]
            τ_k(v) for every vertex v in the graph.
        """

        def extract_subgraph(neighbors: set) -> dict:
            """
            Extracts the neighbor-induced subgraph from the original graph.

            Parameters
            ----------
            neighbors : set
                Set of neighbors of vertex v.

            Returns
            -------
            sub_adj : dict
                Adjacency list representation of the induced subgraph.
            """
            sub_adj = defaultdict(set)
            for u in neighbors:
                for v in G.adj_list[u]:
                    if v in neighbors:
                        sub_adj[u].add(v)
                        sub_adj[v].add(u)
            return sub_adj

        def compute_k_core_components(sub_adj: dict) -> List[set]:
            """
            Given an undirected subgraph, compute its k-core connected components.

            Parameters
            ----------
            sub_adj : dict
                Adjacency list of subgraph.

            Returns
            -------
            components : List[set]
                List of k-core connected components.
            """
            degree = {u: len(sub_adj[u]) for u in sub_adj}
            removed = set()
            queue = deque([u for u in sub_adj if degree[u] < k])

            while queue:
                u = queue.popleft()
                if u in removed:
                    continue
                removed.add(u)
                for v in sub_adj[u]:
                    if v not in removed:
                        degree[v] -= 1
                        if degree[v] == k - 1:
                            queue.append(v)

            def dfs(u, visited, component):
                visited.add(u)
                component.add(u)
                for v in sub_adj[u]:
                    if v not in visited and v not in removed:
                        dfs(v, visited, component)

            visited = set()
            components = []
            for u in sub_adj:
                if u not in visited and u not in removed:
                    comp = set()
                    dfs(u, visited, comp)
                    if comp:
                        components.append(comp)
            return components

        # τ_k(v) 初始化为 0
        τ = [0] * G.vertex_num
        for v in range(G.vertex_num):
            neighbors = set(G.adj_list[v])
            if not neighbors:
                continue
            sub_adj = extract_subgraph(neighbors)
            k_core_components = compute_k_core_components(sub_adj)
            τ[v] = len(k_core_components)

        return τ

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
