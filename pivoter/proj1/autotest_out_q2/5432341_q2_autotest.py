#!/usr/bin/env python3
# Auto-generated for 5432341

STUDENT_ID = "5432341"
STUDENT_NAME = "Yanbo Wang"

# ======= 学生代码 =======
from collections import deque

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
        core = kCoreBaseStructuralDiversity._compute_core_number(G.adj_list, n)

        for v in range(n):
            neighbors = set(G.adj_list[v])
            if not neighbors:
                continue

            if all(core[u] < k for u in neighbors):
                continue

            nbr_graph = {}
            for u in neighbors:
                nbr_graph[u] = []
                for w in G.adj_list[u]:
                    if w in neighbors:
                        nbr_graph[u].append(w)

            k_core_nodes = kCoreBaseStructuralDiversity._compute_k_core(nbr_graph, k)
            if not k_core_nodes:
                continue

            comps = kCoreBaseStructuralDiversity._find_components(nbr_graph, k_core_nodes)
            sd[v] = len(comps)

        return sd

    @staticmethod
    def _compute_core_number(adj_list, n):
        degree = [len(adj_list[i]) for i in range(n)]
        core = degree[:]
        visited = [False] * n
        queue = deque(sorted(range(n), key=lambda x: core[x]))

        while queue:
            v = queue.popleft()
            if visited[v]:
                continue
            visited[v] = True
            for u in adj_list[v]:
                if not visited[u] and core[u] > core[v]:
                    core[u] -= 1
        return core

    @staticmethod
    def _compute_k_core(graph, k):
        deg = {v: len(graph[v]) for v in graph}
        queue = deque([v for v in graph if deg[v] < k])
        removed = set(queue)
        iters = 0

        while queue:
            iters += 1

            v = queue.popleft()
            for u in graph[v]:
                if u not in removed:
                    deg[u] -= 1
                    if deg[u] < k:
                        queue.append(u)
                        removed.add(u)

        return set(graph.keys()) - removed

    @staticmethod
    def _find_components(graph, valid_nodes):
        visited = set()
        components = []

        for v in valid_nodes:
            if v in visited:
                continue
            comp = set()
            queue = deque([v])
            visited.add(v)

            while queue:
                u = queue.popleft()
                comp.add(u)
                for w in graph[u]:
                    if w in valid_nodes and w not in visited:
                        visited.add(w)
                        queue.append(w)

            components.append(comp)

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
