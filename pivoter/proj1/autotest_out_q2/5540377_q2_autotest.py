#!/usr/bin/env python3
# Auto-generated for 5540377

STUDENT_ID = "5540377"
STUDENT_NAME = "Anbang Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque,defaultdict
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
        for v in range(n):
          sd[v] = kCoreBaseStructuralDiversity._compute_k_core_diversity(G, v, k)
        return sd


    @staticmethod
    def _compute_k_core_diversity(G, vertex, k):
        """
        Compute k-core base structural diversity for a single vertex
        """
        # Get neighbors of the vertex
        neighbors = set(G.adj_list[vertex])

        if len(neighbors) == 0:
            return 0

        # Build adjacency list for neighbor-induced subgraph
        neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbors)}
        neighbor_adj = [[] for _ in range(len(neighbors))]

        # Build edges between neighbors
        for neighbor in neighbors:
            for adj_vertex in G.adj_list[neighbor]:
                if adj_vertex in neighbors and adj_vertex != neighbor:
                    neighbor_adj[neighbor_to_idx[neighbor]].append(neighbor_to_idx[adj_vertex])

        # Find all k-cores in the neighbor-induced subgraph
        return kCoreBaseStructuralDiversity._count_k_cores(neighbor_adj, k)

    @staticmethod
    def _count_k_cores(adj_list, k):
        """
        Count the number of k-cores in the graph
        """
        n = len(adj_list)
        if n == 0:
            return 0

        # Calculate degree for each vertex
        degrees = [len(adj_list[i]) for i in range(n)]

        # Use queue for k-core decomposition
        queue = deque()
        removed = [False] * n

        # Add vertices with degree less than k to queue
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)
                removed[i] = True

        # Remove vertices with degree less than k
        while queue:
            v = queue.popleft()
            for neighbor in adj_list[v]:
                if not removed[neighbor]:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        queue.append(neighbor)
                        removed[neighbor] = True

        # Find remaining vertices (vertices in k-core)
        remaining_vertices = [i for i in range(n) if not removed[i]]

        if not remaining_vertices:
            return 0

        # Count the number of connected components
        return kCoreBaseStructuralDiversity._count_connected_components(adj_list, remaining_vertices)

    @staticmethod
    def _count_connected_components(adj_list, vertices):
        """
        Count the number of connected components in the given vertex set
        """
        if not vertices:
            return 0

        vertex_set = set(vertices)
        visited = set()
        components = 0

        for vertex in vertices:
            if vertex not in visited:
                # BFS traversal for connected component
                queue = deque([vertex])
                visited.add(vertex)

                while queue:
                    current = queue.popleft()
                    for neighbor in adj_list[current]:
                        if neighbor in vertex_set and neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

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
