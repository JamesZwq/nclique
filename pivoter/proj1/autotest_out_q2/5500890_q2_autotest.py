#!/usr/bin/env python3
# Auto-generated for 5500890

STUDENT_ID = "5500890"
STUDENT_NAME = "Tao Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core based structural diversity for all vertices in graph G.

        For each vertex v:
        1. Get neighbor set N(v)
        2. Build induced subgraph from N(v)
        3. Find k-cores in the induced subgraph
        4. Return the count of k-cores as τ_k(v)

        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        n = G.vertex_num
        result = []

        for v in range(n):
            # Get neighbors of vertex v
            neighbors = list(G.adj_list[v])

            if len(neighbors) == 0:
                result.append(0)
                continue

            # Special case: if k > max possible degree in neighbor subgraph
            if k > len(neighbors) - 1:
                result.append(0)
                continue

            # Build neighbor-induced subgraph
            neighbor_set = set(neighbors)
            sub_n = len(neighbors)

            # Build adjacency matrix for induced subgraph for easier manipulation
            sub_adj_matrix = [[False] * sub_n for _ in range(sub_n)]
            neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbors)}

            # Add edges between neighbors
            for i, u in enumerate(neighbors):
                for neighbor_of_u in G.adj_list[u]:
                    if neighbor_of_u in neighbor_set and neighbor_of_u != u:
                        j = neighbor_to_idx[neighbor_of_u]
                        sub_adj_matrix[i][j] = True
                        sub_adj_matrix[j][i] = True  # Undirected graph

            # Convert to adjacency list
            sub_adj = [[] for _ in range(sub_n)]
            for i in range(sub_n):
                for j in range(sub_n):
                    if sub_adj_matrix[i][j]:
                        sub_adj[i].append(j)

            # Find k-cores in the induced subgraph
            k_cores_count = kCoreBaseStructuralDiversity._count_k_cores(sub_adj, k)
            result.append(k_cores_count)

        return result

    @staticmethod
    def _count_k_cores(adj_list, k):
        """
        Count the number of k-cores in the given graph.

        A k-core is a maximal connected subgraph where each vertex has degree >= k.

        Parameters
        ----------
        adj_list : List[List[int]]
            Adjacency list representation of the graph
        k : int
            The k value for k-core

        Returns
        -------
        int
            Number of k-cores (connected components where each vertex has degree >= k)
        """
        n = len(adj_list)
        if n == 0:
            return 0

        # Special case: k = 0, all vertices form k-cores, count connected components
        if k == 0:
            visited = [False] * n
            components = 0

            for i in range(n):
                if not visited[i]:
                    components += 1
                    # BFS to mark all connected vertices
                    queue = deque([i])
                    visited[i] = True

                    while queue:
                        curr = queue.popleft()
                        for neighbor in adj_list[curr]:
                            if not visited[neighbor]:
                                visited[neighbor] = True
                                queue.append(neighbor)

            return components

        # For k > 0, apply k-core decomposition
        # Initialize active vertices and degrees
        active = [True] * n
        degrees = [len(adj_list[i]) for i in range(n)]

        # Remove vertices with degree < k iteratively
        changed = True
        while changed:
            changed = False
            for i in range(n):
                if active[i] and degrees[i] < k:
                    active[i] = False
                    changed = True

                    # Update degrees of active neighbors
                    for neighbor in adj_list[i]:
                        if active[neighbor]:
                            degrees[neighbor] -= 1

        # Count connected components in remaining active vertices
        visited = [False] * n
        components = 0

        for i in range(n):
            if active[i] and not visited[i]:
                components += 1
                # BFS to mark all connected active vertices
                queue = deque([i])
                visited[i] = True

                while queue:
                    curr = queue.popleft()
                    for neighbor in adj_list[curr]:
                        if active[neighbor] and not visited[neighbor]:
                            visited[neighbor] = True
                            queue.append(neighbor)

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
