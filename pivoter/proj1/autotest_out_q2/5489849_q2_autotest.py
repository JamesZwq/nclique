#!/usr/bin/env python3
# Auto-generated for 5489849

STUDENT_ID = "5489849"
STUDENT_NAME = "Ating Chen"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules.
################################################################################
from collections import deque, defaultdict

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Calculate the k-core-based structural diversity for each vertex in the graph.

        Input:
          G: UndirectedUnweightedGraph object with adj_list attribute
             Or dict, undirected unweighted graph in adjacency list, G[u] = [v1, v2, ...]
          k: int, k-core parameter

        Output:
          An array tau of length |V|, where tau[v] is the k-core structural diversity of vertex v
        """
        # Handle both graph representations (UndirectedUnweightedGraph or dict)
        if hasattr(G, 'adj_list'):
            adj_list = G.adj_list
            n = G.vertex_num
        else:
            adj_list = G
            n = max(G.keys()) + 1 if G else 0

        # Initialize result array
        tau = [0] * n

        # Compute k-core vertices
        k_core_vertices = kCoreBaseStructuralDiversity._compute_k_core(adj_list, k)

        # If k-core is empty, return all zeros
        if not k_core_vertices:
            return tau

        # Pre-compute a set version of the graph for faster lookups
        adj_sets = {}
        for v in range(n) if hasattr(G, 'adj_list') else G:
            if v < len(adj_list):
                adj_sets[v] = set(adj_list[v])

        # For each vertex, compute its structural diversity
        for v in range(n) if hasattr(G, 'adj_list') else G:
            if v in k_core_vertices:
                # Get neighbors that are also in the k-core
                k_core_neighbors = set(u for u in adj_list[v] if u in k_core_vertices)

                # If vertex has enough neighbors in the k-core
                if len(k_core_neighbors) >= k:
                    # Count connected components in the induced subgraph
                    tau[v] = kCoreBaseStructuralDiversity._count_components(
                        adj_sets, k_core_neighbors, k, k_core_vertices)

        return tau

    @staticmethod
    def _compute_k_core(adj_list, k):
        """
        Compute the k-core of the graph.

        Returns a set of vertices in the k-core.
        """
        # Handle both list and dict graph representations
        if isinstance(adj_list, list):
            vertices = range(len(adj_list))
        else:
            vertices = adj_list.keys()

        # Initialize degrees and neighbors
        degrees = {}
        neighbors = {}

        for v in vertices:
            if v < len(adj_list) if isinstance(adj_list, list) else True:
                neighbors[v] = set(adj_list[v])
                degrees[v] = len(neighbors[v])

        # Queue vertices with degree < k for removal
        queue = deque([v for v in degrees if degrees[v] < k])
        removed = set()

        while queue:
            v = queue.popleft()
            removed.add(v)

            # Process each neighbor of v
            for u in list(neighbors[v]):
                if u not in removed and u in neighbors:
                    # Remove edge (u, v)
                    neighbors[u].discard(v)
                    # Update degree
                    degrees[u] -= 1
                    # If degree drops below k, queue for removal
                    if degrees[u] < k:
                        queue.append(u)

        # Return vertices that were not removed
        k_core_vertices = set(v for v in degrees if v not in removed)

        # Verify all remaining vertices have degree >= k
        for v in list(k_core_vertices):
            valid_neighbors = [u for u in neighbors[v] if u in k_core_vertices]
            if len(valid_neighbors) < k:
                k_core_vertices.remove(v)

        return k_core_vertices

    @staticmethod
    def _count_components(adj_sets, vertices, k, k_core_vertices):
        """
        Count connected components in the k-core of the induced subgraph on vertices.

        Args:
            adj_sets: Graph with adjacency lists as sets for faster lookups
            vertices: Set of vertices to consider (neighbors of a vertex)
            k: k-core parameter
            k_core_vertices: Pre-computed set of vertices in the k-core
        """
        if len(vertices) < k:
            return 0

        # Create the induced subgraph of vertices that are also in the k-core
        vertices = vertices.intersection(k_core_vertices)
        if len(vertices) < k:
            return 0

        # Build the induced subgraph
        subgraph = {}
        for u in vertices:
            subgraph[u] = vertices.intersection(adj_sets[u])

        # Apply k-core decomposition to the induced subgraph
        degrees = {u: len(subgraph[u]) for u in subgraph}
        to_remove = [u for u in subgraph if degrees[u] < k]

        while to_remove:
            # Remove vertices with degree < k
            for u in to_remove:
                for v in subgraph[u]:
                    subgraph[v].remove(u)
                    degrees[v] -= 1
                del subgraph[u]
                del degrees[u]

            # Find new vertices to remove
            to_remove = [u for u in subgraph if degrees[u] < k]

        # Count connected components using BFS
        if not subgraph:
            return 0

        visited = set()
        components = 0

        for u in subgraph:
            if u not in visited:
                components += 1
                queue = deque([u])
                visited.add(u)

                while queue:
                    node = queue.popleft()
                    for neighbor in subgraph[node]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

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
