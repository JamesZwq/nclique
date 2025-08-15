#!/usr/bin/env python3
# Auto-generated for 5458019

STUDENT_ID = "5458019"
STUDENT_NAME = "Jung-An Lin"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################

from collections import deque, defaultdict


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

        def extract_neighbor_induced_subgraph(graph, N_v):
            subgraph = defaultdict(set)
            for u in N_v:
                subgraph[u]
                for v in graph[u]:
                    if v in N_v:
                        subgraph[u].add(v)
            return subgraph


        def perform_k_core_decomposition(graph, k):
            def deg(v):
                return len(graph[v])

            changed = True
            while changed:
                changed = False
                for v in list(graph.keys()):
                    if deg(v) < k:
                        for nb in graph[v]:
                            graph[nb].discard(v)
                            changed = True
                        graph.pop(v)
            return graph


        def identify_connected_components(graph):
            visited = set()
            components = []

            def bfs(s):
                queue = deque([s])
                visited.add(s)
                component = []

                while queue:
                    v = queue.popleft()
                    component.append(v)
                    for nb in graph[v]:
                        if nb not in visited:
                            queue.append(nb)
                            visited.add(nb)
                return component

            for v in graph:
                if v not in visited:
                    component = bfs(v)
                    components.append(component)
            return components


        def compute_tau_k_v(graph, k, v):
            N_v = set(graph[v])
            nbr_v = extract_neighbor_induced_subgraph(graph, N_v)
            if k == 0:
                core_k = nbr_v
            else:
                core_k = perform_k_core_decomposition(nbr_v, k)
            components = identify_connected_components(core_k)
            return len(components)


        n = G.vertex_num
        graph = G.adj_list
        T = list()
        for v in range(n):
            tau = compute_tau_k_v(graph, k, v)
            T.append(tau)
        return T


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
