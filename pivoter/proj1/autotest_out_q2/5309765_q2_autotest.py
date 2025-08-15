#!/usr/bin/env python3
# Auto-generated for 5309765

STUDENT_ID = "5309765"
STUDENT_NAME = "Dharshini Loganathan"

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

        n = G.vertex_num  # total number of vertices in the graph
        tau_k = [0] * n   # initialize tau_k list with zeros

        # Function to get the subgraph induced by the neighbors of vertex v
        def get_neighbour_subgraph(v):
            neighbors = set(G.adj_list[v])  # get all neighbors of vertex v
            sub_edges = []

            # loop over each neighbor and collect edges between neighbors
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbors:
                        sub_edges.append((u, w))

            return neighbors, sub_edges  # return nodes and edges of subgraph

        # Function to compute how many connected components exist in k-core
        def compute_k_core_components(nodes, edges):
            # build adjacency list for the subgraph
            adj = {u: set() for u in nodes}
            for u, v in edges:
                adj[u].add(v)
                adj[v].add(u)

            # calculate degrees of each node
            deg = {u: len(adj[u]) for u in nodes}
            # put all nodes with degree < k into a queue
            q = deque([u for u in nodes if deg[u] < k])
            alive = set(nodes)  # nodes still in the k-core

            # remove all nodes that can't stay in the k-core
            while q:
                u = q.popleft()
                if u not in alive:
                    continue
                alive.remove(u)
                for v in adj[u]:
                    if v in alive:
                        deg[v] -= 1
                        if deg[v] < k:
                            q.append(v)

            # count number of connected components in remaining nodes
            visited = set()
            components = 0

            for u in alive:
                if u not in visited:
                    components += 1
                    stack = [u]
                    while stack:
                        node = stack.pop()
                        if node in visited:
                            continue
                        visited.add(node)
                        for nei in adj[node]:
                            if nei in alive and nei not in visited:
                                stack.append(nei)

            return components

        # For every vertex in the graph, compute tau_k value
        for v in range(n):
            nodes, edges = get_neighbour_subgraph(v)  # get neighbor subgraph
            if nodes:
                tau_k[v] = compute_k_core_components(nodes, edges)

        return tau_k  # return final list of tau_k values for each vertex

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
