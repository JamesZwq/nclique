#!/usr/bin/env python3
# Auto-generated for 5521365

STUDENT_ID = "5521365"
STUDENT_NAME = "Zixi Hong"

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
        sd = [0] * n  # structural diversity array

        for v in range(n):
            neighbors = G.adj_list[v]  # directly access adjacency list

            # build neighbor-induced subgraph (adjacency list)
            sub_adj = {}
            for u in neighbors:
                sub_adj[u] = [w for w in G.adj_list[u] if w in neighbors]

            # Calculate the number of connected components in the k-core of this subgraph
            sd[v] = kCoreBaseStructuralDiversity.compute_k_core_components(sub_adj, k)

        return sd

    @staticmethod
    def compute_k_core_components(adj, k):
        """
        adj: dict[int, list[int]] adjacency list (subgraph)
        k: int k-core threshold
        return: number of connected components in remaining k-core graph
        """
        degree = {u: len(adj[u]) for u in adj}
        deleted = set()
        queue = deque([u for u in adj if degree[u] < k])

        while queue:
            u = queue.popleft()
            if u in deleted:
                continue
            deleted.add(u)

            for v in adj[u]:
                if v not in deleted:
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)

        # DFS statistics connected components
        visited = set()
        components = 0

        def dfs(start):
            stack = [start]
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                for neighbor in adj[node]:
                    if neighbor not in visited and neighbor not in deleted:
                        stack.append(neighbor)

        for u in adj:
            if u not in deleted and u not in visited:
                dfs(u)
                components += 1

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
