#!/usr/bin/env python3
# Auto-generated for 5570806

STUDENT_ID = "5570806"
STUDENT_NAME = "Yuanyuan Tang"

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

        # Calculate τ_k for each node separately
        for node in range(n):
            neighbors = G.adj_list[node]
            ## Special case 1: The current node has no neighbors
            if not neighbors:
                sd[node] = 0
                continue

            ## Special case 2: k == 0
            if k == 0:
                sub_adj = kCoreBaseStructuralDiversity._build_subgraph(G, neighbors)
                sd[node] = kCoreBaseStructuralDiversity._count_components(sub_adj, sub_adj.keys())
                continue

            ## Special case 3: The number of neighbors is less than k
            if len(neighbors) < k:
                sd[node] = 0
                continue

            # 1. Construct a neighbor induction subgraph
            sub_adj = kCoreBaseStructuralDiversity._build_subgraph(G, neighbors)

            # 2. Peeling treatment
            remaining = kCoreBaseStructuralDiversity._peeling(sub_adj, k)

            ## Special case 4: There are no remaining nodes after peeling
            if not remaining:
                sd[node] = 0
                continue

            # 3. Count the number of connected components formed by the remaining nodes in k-core
            components = kCoreBaseStructuralDiversity._count_components(sub_adj, remaining)
            sd[node] = components

        return sd

    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################

    @staticmethod
    # Construct a neighbor induction subgraph
    def _build_subgraph(G, neighbors):
        S = set(neighbors)
        sub_adj = {u: [] for u in S}
        for u in S:
            for v in G.adj_list[u]:
                if v in S:
                    sub_adj[u].append(v)
        return sub_adj

    @staticmethod
    # Execute the peeling logic and return the remaining nodes
    def _peeling(sub_adj, k):
        degree = {u: len(sub_adj[u]) for u in sub_adj}
        removed = set()
        queue = deque()

        # Initially, add all nodes with a degree less than k to the queue
        for u in sub_adj:
            if degree[u] < k:
                queue.append(u)
                removed.add(u)

        # BFS-like
        while queue:
            current = queue.popleft()
            for nb in sub_adj[current]:
                if nb not in removed:
                    degree[nb] -= 1
                    if degree[nb] < k:
                        queue.append(nb)
                        removed.add(nb)

        remaining = [u for u in sub_adj if u not in removed]
        return remaining

    @staticmethod
    # Statistical connected component
    def _count_components(sub_adj, remaining_nodes):
        visited = set()
        components = 0
        queue = deque()

        for u in remaining_nodes:
            if u not in visited:
                components += 1
                queue.append(u)
                visited.add(u)
                while queue:
                    current = queue.popleft()
                    for nb in sub_adj[current]:
                        if nb not in visited and nb in remaining_nodes:
                            visited.add(nb)
                            queue.append(nb)

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
