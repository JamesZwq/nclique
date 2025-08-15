#!/usr/bin/env python3
# Auto-generated for 5541704

STUDENT_ID = "5541704"
STUDENT_NAME = "zheng Guan"

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

        # Process each vertex
        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                # If no neighbors, τ_k(v) = 0
                sd[v] = 0
                continue

            # Count the number of k-cores in the subgraph induced by the neighbors of v
            sd[v] = kCoreBaseStructuralDiversity._count_k_cores(G, neighbors, k)

        return sd

    @staticmethod
    def _count_k_cores(G, vertices, k):
        """
        Count the number of k-cores in the subgraph induced by a given set of vertices.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The original graph.
        vertices : List[int]
            The set of vertices to induce the subgraph.
        k : int
            The core level to compute.

        Returns
        -------
        int
            Number of k-cores (connected components in the k-core subgraph).
        """
        if k == 0:
            # For k=0, simply return the number of connected components
            return kCoreBaseStructuralDiversity._count_connected_components(G, vertices)

        # Build induced subgraph
        vertex_set = set(vertices)
        vertex_to_idx = {v: i for i, v in enumerate(vertices)}
        n = len(vertices)

        adj = [[] for _ in range(n)]
        degree = [0] * n

        # Build adjacency list and degree list for the induced subgraph
        for i, u in enumerate(vertices):
            for v in G.adj_list[u]:
                if v in vertex_set:
                    j = vertex_to_idx[v]
                    adj[i].append(j)
                    degree[i] += 1

        # Perform k-core decomposition by iteratively removing low-degree nodes
        removed = [False] * n
        queue = deque()

        for i in range(n):
            if degree[i] < k:
                queue.append(i)
                removed[i] = True

        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if not removed[v]:
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)
                        removed[v] = True

        # Remaining vertices after decomposition
        remaining = [i for i in range(n) if not removed[i]]
        if not remaining:
            return 0

        # Count the number of connected components in the remaining graph
        visited = [False] * n
        k_core_count = 0

        for start in remaining:
            if not visited[start]:
                # BFS to explore a connected component
                queue = deque([start])
                visited[start] = True

                while queue:
                    u = queue.popleft()
                    for v in adj[u]:
                        if not removed[v] and not visited[v]:
                            visited[v] = True
                            queue.append(v)

                k_core_count += 1

        return k_core_count

    @staticmethod
    def _count_connected_components(G, vertices):
        """
        Count the number of connected components in the subgraph induced by a set of vertices.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The original graph.
        vertices : List[int]
            The set of vertices to induce the subgraph.

        Returns
        -------
        int
            Number of connected components in the induced subgraph.
        """
        vertex_set = set(vertices)
        visited = set()
        component_count = 0

        for start in vertices:
            if start not in visited:
                # BFS to explore a connected component
                queue = deque([start])
                visited.add(start)

                while queue:
                    u = queue.popleft()
                    for v in G.adj_list[u]:
                        if v in vertex_set and v not in visited:
                            visited.add(v)
                            queue.append(v)

                component_count += 1

        return component_count

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
