#!/usr/bin/env python3
# Auto-generated for 5360671

STUDENT_ID = "5360671"
STUDENT_NAME = "Guandong Wang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

# Please see Q2.pdf for a more detailed explanation for my algorithm :)
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # This function constructs the neighbour induced subgraph from a given set of nodes
    @staticmethod
    def build_induced_subgraph(nodes, G):
        node_set = set(nodes)
        subgraph = {u: set() for u in nodes}
        for u in nodes:
            for v in G.adj_list[u]:
                if v in node_set:
                    subgraph[u].add(v)
        return subgraph

    # Function to find all connected components in an undirected graph using BFS
    @staticmethod
    def get_connected_components(graph):
        visited = set()
        components = []

        for node in graph:

            if node in visited:
                continue

            queue = deque([node])
            visited.add(node)
            comp = []

            while queue:
                u = queue.popleft()
                comp.append(u)
                for v in graph[u]:
                    if v not in visited:
                        visited.add(v)
                        queue.append(v)
            components.append(comp)

        return components



    # This function is the algorithm for k-core similar to the lecture example. It peels a graph to extract its k-core subgraph using a degree based algorithm.
    @staticmethod
    def peel_kcore(graph, k):
        degrees = {u: len(graph[u]) for u in graph}
        remaining = set(graph)
        queue = deque([u for u in graph if degrees[u] < k])

        while queue:
            u = queue.popleft()
            remaining.discard(u)
            for v in graph[u]:
                if v in remaining:
                    degrees[v] -= 1
                    if degrees[v] < k and v not in queue:
                        queue.append(v)
        return remaining


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

        for v in range(n):
            nbrs = G.adj_list[v]

            # Edge case checker for vertex with no neighbours
            if not nbrs:
                continue

            # Build the induced subgraphs from the neighbours of v
            subgraph = kCoreBaseStructuralDiversity.build_induced_subgraph(nbrs, G)

            # Apply k-core peeling to extract the remaining nodes
            remaining = kCoreBaseStructuralDiversity.peel_kcore(subgraph, k)

            # Rebuild subgraph using remaining nodes and count connected components
            if remaining:
                pruned_graph = kCoreBaseStructuralDiversity.build_induced_subgraph(remaining, G)
                components = kCoreBaseStructuralDiversity.get_connected_components(pruned_graph)
                sd[v] = len(components)

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
