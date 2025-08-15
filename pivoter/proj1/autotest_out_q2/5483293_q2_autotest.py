#!/usr/bin/env python3
# Auto-generated for 5483293

STUDENT_ID = "5483293"
STUDENT_NAME = "Yanko Zhang"

# ======= 学生代码 =======
################################################################################
# This algorithm first extracts the neighbor-induced subgraph of each node v 
# in the graph, and extracts the k-core structure in the subgraph by iteratively 
# stripping off nodes with a degree less than k. Then, the k-core subgraph 
# obtained by BFS is traversed, and the number of connected components in it is 
# counted as the structural diversity value of node v.
# Time complexity analysis:
# Let
# Number of vertices n
# Number of edges m
# The degree of each vertex d
# The max degree of each vertex D
# The algorithm constructs a neighbor induced subgraph for each node v, extracts 
# k-cores from it and calculates the number of connected components. 
# For a single node, the time complexity of constructing the neighbor subgraph is 
# O(dD), and the total cost of extracting k-cores and traversing connected 
# components is O(d^2), where d is the degree of node v and D is the maximum degree 
# in the graph. 
# Therefore, the time complexity of the overall algorithm is O(Dm+∑d^2).
################################################################################
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
        tau = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if len(neighbors) == 0:
                continue

            # Build subgraph
            subgraph_nodes = set(neighbors)
            subgraph_adj = {u: [] for u in subgraph_nodes}

            for u in subgraph_nodes:
                for w in G.adj_list[u]:
                    if w in subgraph_nodes:
                        subgraph_adj[u].append(w)

            # Remove nodes with degree < k
            def extract_k_core(adj_dict):
                deg = {u: len(adj_dict[u]) for u in adj_dict}
                queue = deque([u for u in deg if deg[u] < k])
                while queue:
                    u = queue.popleft()
                    for v in adj_dict[u]:
                        if v in deg:
                            deg[v] -= 1
                            if deg[v] == k - 1:
                                queue.append(v)
                    deg.pop(u)
                # Build k-core subgraph
                core_adj = {}
                for u in deg:
                    core_adj[u] = [v for v in adj_dict[u] if v in deg]
                return core_adj

            k_core_adj = extract_k_core(subgraph_adj)

            # Count connected components in k-core subgraph
            visited = set()
            def bfs(start):
                queue = deque([start])
                visited.add(start)
                while queue:
                    u = queue.popleft()
                    for v in k_core_adj[u]:
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)

            count = 0
            for u in k_core_adj:
                if u not in visited:
                    bfs(u)
                    count += 1
            tau[v] = count

        return tau

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
