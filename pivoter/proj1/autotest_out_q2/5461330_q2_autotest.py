#!/usr/bin/env python3
# Auto-generated for 5461330

STUDENT_ID = "5461330"
STUDENT_NAME = "Nora Hong"

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

        # For each vertex v, compute τ_k(v)
        for v in range(n):
            # Get neighbors of v
            neighbors = G.adj_list[v]

            if len(neighbors) == 0:
                tau[v] = 0
                continue

            # Build neighbor-induced subgraph
            # Create mapping from original vertex id to index in neighbor list
            neighbor_to_idx = {}
            for i, neighbor in enumerate(neighbors):
                neighbor_to_idx[neighbor] = i

            m = len(neighbors)

            # Build adjacency list for the neighbor-induced subgraph
            nbr_adj = [[] for _ in range(m)]
            for i, u in enumerate(neighbors):
                for w in G.adj_list[u]:
                    if w in neighbor_to_idx and w != u:  # w is also a neighbor of v, and not u itself
                        j = neighbor_to_idx[w]
                        nbr_adj[i].append(j)

            # Find number of k-cores in the neighbor-induced subgraph
            tau[v] = kCoreBaseStructuralDiversity._count_k_cores(nbr_adj, k)

        return tau

    @staticmethod
    def _count_k_cores(adj_list, k):
        """
        Count the number of k-cores (connected components that form k-cores)
        in a graph represented by adjacency list.
        """
        n = len(adj_list)
        if n == 0:
            return 0

        # Special case for k=0: every connected component forms a 0-core
        if k == 0:
            return kCoreBaseStructuralDiversity._count_connected_components(adj_list)

        # For k > 0, find k-core using iterative vertex removal
        degrees = [len(adj_list[i]) for i in range(n)]
        removed = [False] * n
        queue = deque()

        # Initially add all vertices with degree < k to queue
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)

        # Iteratively remove vertices with degree < k
        while queue:
            u = queue.popleft()
            if removed[u]:
                continue
            removed[u] = True

            # Update degrees of neighbors
            for v in adj_list[u]:
                if not removed[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)

        # Build adjacency list for remaining vertices (k-core vertices)
        remaining_vertices = [i for i in range(n) if not removed[i]]
        if not remaining_vertices:
            return 0

        # Create mapping from original index to new index
        old_to_new = {}
        for i, v in enumerate(remaining_vertices):
            old_to_new[v] = i

        # Build new adjacency list with only remaining vertices
        k_core_adj = [[] for _ in range(len(remaining_vertices))]
        for i, u in enumerate(remaining_vertices):
            for v in adj_list[u]:
                if not removed[v]:  # v is also in k-core
                    k_core_adj[i].append(old_to_new[v])

        # Count connected components in k-core
        return kCoreBaseStructuralDiversity._count_connected_components(k_core_adj)

    @staticmethod
    def _count_connected_components(adj_list):
        """Count the number of connected components in a graph."""
        n = len(adj_list)
        if n == 0:
            return 0

        visited = [False] * n
        components = 0

        for i in range(n):
            if not visited[i]:
                components += 1
                # BFS to mark all vertices in this component
                queue = deque([i])
                visited[i] = True

                while queue:
                    u = queue.popleft()
                    for v in adj_list[u]:
                        if not visited[v]:
                            visited[v] = True
                            queue.append(v)

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
