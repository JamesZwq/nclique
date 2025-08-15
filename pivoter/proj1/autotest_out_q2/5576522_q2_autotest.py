#!/usr/bin/env python3
# Auto-generated for 5576522

STUDENT_ID = "5576522"
STUDENT_NAME = "Shuli Zeng"

# ======= 学生代码 =======
from collections import deque, defaultdict
import time
import os


class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core-based structural diversity for all vertices.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The input graph with G.vertex_num vertices and G.adj_list adjacency list
        k : int
            The k value for k-core computation

        Returns
        -------
        List[int]
            Array of τ_k(v) values for all vertices v
        """
        n = G.vertex_num
        tau = [0] * n

        # Handle edge cases
        if k < 0 or n == 0:
            return tau

        # Process each vertex
        for v in range(n):
            # Get neighbors of v
            neighbors = list(G.adj_list[v])

            # If not enough neighbors for k-core, τ_k(v) = 0
            if len(neighbors) < k:
                tau[v] = 0
                continue

            # Compute τ_k(v) = number of k-cores in neighbor-induced subgraph
            tau[v] = kCoreBaseStructuralDiversity._count_k_cores_in_subgraph(
                G.adj_list, neighbors, k
            )

        return tau

    @staticmethod
    def _count_k_cores_in_subgraph(adj_list, vertices, k):
        """
        Count the number of k-cores in the subgraph induced by given vertices.

        Parameters
        ----------
        adj_list : adjacency list of the original graph
        vertices : list of vertices that induce the subgraph
        k : k value for k-core

        Returns
        -------
        int : number of k-cores (connected components in k-core)
        """
        if len(vertices) < k:
            return 0

        # Create vertex set for quick lookup
        vertex_set = set(vertices)

        # Build degree map for vertices in the induced subgraph
        degree = {}
        for v in vertices:
            # Count neighbors that are also in the vertex set
            deg = 0
            for u in adj_list[v]:
                if u in vertex_set:
                    deg += 1
            degree[v] = deg

        # Find k-core vertices using peeling algorithm
        removed = set()
        queue = deque()

        # Initialize queue with vertices having degree < k
        for v in vertices:
            if degree[v] < k:
                queue.append(v)
                removed.add(v)

        # Peeling process
        while queue:
            v = queue.popleft()

            # Update degrees of neighbors
            for u in adj_list[v]:
                if u in vertex_set and u not in removed:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        removed.add(u)

        # Remaining vertices form the k-core
        k_core_vertices = [v for v in vertices if v not in removed]

        if not k_core_vertices:
            return 0

        # Count connected components in k-core
        return kCoreBaseStructuralDiversity._count_connected_components(
            adj_list, k_core_vertices
        )

    @staticmethod
    def _count_connected_components(adj_list, vertices):
        """
        Count connected components in the subgraph induced by given vertices.

        Parameters
        ----------
        adj_list : adjacency list of the original graph
        vertices : list of vertices that induce the subgraph

        Returns
        -------
        int : number of connected components
        """
        if not vertices:
            return 0

        vertex_set = set(vertices)
        visited = set()
        component_count = 0

        # Use DFS to find connected components
        for start in vertices:
            if start not in visited:
                component_count += 1

                # DFS from start vertex
                stack = [start]
                visited.add(start)

                while stack:
                    v = stack.pop()

                    # Visit all unvisited neighbors in the vertex set
                    for u in adj_list[v]:
                        if u in vertex_set and u not in visited:
                            visited.add(u)
                            stack.append(u)

        return component_count

# ==================== Graph Class Definition ====================
class UndirectedUnweightedGraph:
    """Graph class that matches the expected interface"""
    def __init__(self, vertex_num):
        self.vertex_num = vertex_num
        self.adj_list = [[] for _ in range(vertex_num)]

    def add_edge(self, u, v):
        """Add an undirected edge between u and v"""
        self.adj_list[u].append(v)
        self.adj_list[v].append(u)

# ==================== Helper Functions ====================
# def load_graph(filename):
#     """Load graph from file"""
#     print(f" Loading graph: {filename}")

#     with open(filename, 'r') as f:
#         # Read first line: num_vertices num_edges
#         first_line = f.readline().strip().split()
#         num_vertices = int(first_line[0])
#         num_edges = int(first_line[1])

#         # Create graph
#         G = UndirectedUnweightedGraph(num_vertices)

#         # Read edges
#         edges_read = 0
#         for line in f:
#             if line.strip():
#                 parts = line.strip().split()
#                 if len(parts) >= 2:
#                     u, v = int(parts[0]), int(parts[1])
#                     G.add_edge(u, v)
#                     edges_read += 1

#         print(f"   ✓ Loaded: {num_vertices} vertices, {edges_read} edges")

#         return G

# def load_answer(filename):
#     """Load expected answer from file"""
#     with open(filename, 'r') as f:
#         content = f.read().strip()

#     # Parse answer values
#     values = list(map(int, content.split()))

#     return values

# # ==================== Verification Function ====================
# def verify_all_tests():
#     """Verify all test cases with corrected comparison"""
#     print(" Q2 k-Core Based Structural Diversity - Verification")
#     print("=" * 70)

#     test_cases = [
#         ('/content/dataset_9312/data_15.graph.txt', 6, '/content/dataset_9312/ans_15_6.txt'),
#         ('/content/dataset_9312/data_15.graph.txt', 9, '/content/dataset_9312/ans_15_9.txt'),
#         ('/content/dataset_9312/data_20.graph.txt', 6, '/content/dataset_9312/ans_20_6.txt'),
#         ('/content/dataset_9312/data_20.graph.txt', 7, '/content/dataset_9312/ans_20_7.txt'),
#         ('/content/dataset_9312/data_80.graph.txt', 9, '/content/dataset_9312/ans_80_9.txt'),
#         ('/content/dataset_9312/data_80.graph.txt', 14, '/content/dataset_9312/ans_80_14.txt'),
#     ]

#     passed = 0
#     total = 0

#     for graph_file, k, answer_file in test_cases:
#         if os.path.exists(graph_file) and os.path.exists(answer_file):
#             total += 1

#             # Load and run
#             G = load_graph(graph_file)
#             result = kCoreBaseStructuralDiversity.process(G, k)
#             expected = load_answer(answer_file)

#             # Build answer array
#             tau_dist = defaultdict(int)
#             for tau in result:
#                 tau_dist[tau] += 1

#             max_tau = max(result) if result else 0
#             our_answer = [k]
#             for tau_val in range(max_tau + 1):
#                 our_answer.append(tau_dist[tau_val])

#             # Trim to expected length
#             our_answer = our_answer[:len(expected)]

#             # Check match
#             if our_answer == expected:
#                 print(f" {graph_file}, k={k}: Passed")
#                 passed += 1
#             else:
#                 print(f" {graph_file}, k={k}: Failed")
#                 print(f"   Expected: {expected}")
#                 print(f"   Got: {our_answer}")

#                 # Show detailed τ distribution
#                 print(f"   τ distribution:")
#                 for tau_val in range(max_tau + 1):
#                     if tau_dist[tau_val] > 0:
#                         print(f"     τ_{k} = {tau_val}: {tau_dist[tau_val]} vertices")

#     print(f"\n Summary: {passed}/{total} tests passed")

#     if passed == total:
#         print(" All tests passed! Your algorithm is correct!")

#     return passed == total

# # ==================== Detailed Test Function ====================
# def test_case_detailed(graph_file, k, answer_file=None):
#     """Detailed single test case analysis"""

#     print(f"\n{'='*60}")
#     print(f" Test: {graph_file}, k={k}")
#     print(f"{'='*60}")

#     # Load graph
#     G = load_graph(graph_file)

#     # Run algorithm
#     print(f"\n  Running algorithm...")
#     start_time = time.time()
#     result = kCoreBaseStructuralDiversity.process(G, k)
#     end_time = time.time()

#     print(f"   ✓ Time taken: {end_time - start_time:.4f} seconds")

#     # Calculate τ distribution
#     tau_distribution = defaultdict(int)
#     for tau in result:
#         tau_distribution[tau] += 1

#     # Sort by τ value
#     max_tau = max(result) if result else 0

#     print(f"\n Result Summary:")
#     print(f"   • Total vertices: {G.vertex_num}")
#     print(f"   • Max τ_{k} value: {max_tau}")

#     # Show distribution
#     print(f"\n τ_{k} value distribution:")
#     for tau_val in range(max_tau + 1):
#         count = tau_distribution[tau_val]
#         if count > 0:
#             print(f"   τ_{k} = {tau_val}: {count} vertices")

#     # Load and compare with expected answer if provided
#     if answer_file and os.path.exists(answer_file):
#         print(f"\n Comparing with expected answer...")
#         expected = load_answer(answer_file)

#         # Build our answer in the expected format
#         our_answer = [k]  # Start with k
#         for tau_val in range(max_tau + 1):
#             our_answer.append(tau_distribution[tau_val])

#         # Trim to match expected length
#         our_answer = our_answer[:len(expected)]

#         print(f"   Expected: {expected}")
#         print(f"   Got: {our_answer}")

#         if our_answer == expected:
#             print(f"    Perfect match!")
#         else:
#             print(f"    Mismatch")

# # ==================== Main ====================
# if __name__ == "__main__":
#     # Run all tests
#     verify_all_tests()

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
