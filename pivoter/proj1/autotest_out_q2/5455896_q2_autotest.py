#!/usr/bin/env python3
# Auto-generated for 5455896

STUDENT_ID = "5455896"
STUDENT_NAME = "Tianjie Xu"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        calculator = kCoreBaseStructuralDiversity()

        for v in range(n):
            neighbors = set(G.adj_list[v])

            if not neighbors:
                continue

            # Build neighbor-induced subgraph
            nbr_graph = {u: [w for w in G.adj_list[u] if w in neighbors] for u in neighbors}

            # Count k-cores
            sd[v] = len(calculator.find_kcores(nbr_graph, k))

        return sd

    def find_kcores(self, graph, k):
        if not graph:
            return []

        k_core_vertices = self.compute_kcore(graph, k)
        return self.find_components(graph, k_core_vertices) if k_core_vertices else []

    def compute_kcore(self, graph, k):
        degree = {v: len(graph[v]) for v in graph}
        removed = set()
        queue = deque(v for v in graph if degree[v] < k)
        removed.update(queue)

        while queue:
            v = queue.popleft()
            for u in graph[v]:
                if u not in removed:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        removed.add(u)

        return set(graph.keys()) - removed

    def find_components(self, graph, vertices):
        if not vertices:
            return []

        visited = set()
        components = []

        for v in vertices:
            if v in visited:
                continue

            component = set()
            queue = deque([v])
            component.add(v)
            visited.add(v)

            while queue:
                u = queue.popleft()
                for w in graph[u]:
                    if w in vertices and w not in visited:
                        queue.append(w)
                        component.add(w)
                        visited.add(w)

            components.append(component)

        return components


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
