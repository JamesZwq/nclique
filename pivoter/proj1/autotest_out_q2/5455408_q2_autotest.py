#!/usr/bin/env python3
# Auto-generated for 5455408

STUDENT_ID = "5455408"
STUDENT_NAME = "Yurong Yang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

################################################################################
# Standard Library imports for graph processing
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core based structural diversity for each vertex in graph G

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            Input undirected unweighted graph
        k : int
            Parameter k for k-core computation

        Returns
        -------
        List[int]
            τ_k(v) values representing structural diversity for all vertices
        """

        # Initialize result array with zeros
        structural_diversity = [0] * G.vertex_num

        # Process each vertex to compute its structural diversity
        for vertex_id in range(G.vertex_num):
            neighbor_list = G.adj_list[vertex_id]

            # Early termination for vertices with insufficient neighbors
            if len(neighbor_list) < k:
                continue

            # Build induced subgraph from vertex neighbors
            subgraph_adjacency = kCoreBaseStructuralDiversity._construct_induced_subgraph(
                G, neighbor_list
            )

            # Calculate k-core components in the induced subgraph
            structural_diversity[vertex_id] = kCoreBaseStructuralDiversity._compute_kcore_components(
                neighbor_list, subgraph_adjacency, k
            )

        return structural_diversity

    @staticmethod
    def _construct_induced_subgraph(graph, vertex_neighbors):
        """
        Build adjacency representation for subgraph induced by given vertices

        Args:
            graph: Original graph structure
            vertex_neighbors: List of vertices forming the induced subgraph

        Returns:
            dict: Adjacency mapping for the induced subgraph
        """
        neighbor_set = set(vertex_neighbors)
        induced_edges = defaultdict(set)

        # Iterate through vertices to find internal connections
        for current_vertex in vertex_neighbors:
            adjacent_vertices = graph.adj_list[current_vertex]

            for adjacent_vertex in adjacent_vertices:
                # Include edge only if both endpoints are in neighbor set
                if (adjacent_vertex in neighbor_set and
                    current_vertex < adjacent_vertex):  # Avoid duplicate edges
                    induced_edges[current_vertex].add(adjacent_vertex)
                    induced_edges[adjacent_vertex].add(current_vertex)

        return induced_edges

    @staticmethod
    def _compute_kcore_components(vertices, edge_mapping, k_threshold):
        """
        Find number of k-core connected components in given subgraph

        Args:
            vertices: List of vertices in subgraph
            edge_mapping: Adjacency representation of subgraph
            k_threshold: Minimum degree requirement for k-core

        Returns:
            int: Number of connected k-core components
        """
        if not vertices:
            return 0

        # Apply k-core decomposition to filter vertices
        remaining_vertices = kCoreBaseStructuralDiversity._apply_kcore_filtering(
            vertices, edge_mapping, k_threshold
        )

        # Count connected components in remaining k-core vertices
        return kCoreBaseStructuralDiversity._count_components(
            remaining_vertices, edge_mapping
        )

    @staticmethod
    def _apply_kcore_filtering(vertex_list, adjacency_map, min_degree):
        """
        Remove vertices with degree less than k using iterative approach

        Args:
            vertex_list: Initial set of vertices
            adjacency_map: Edge connections between vertices
            min_degree: Minimum degree threshold

        Returns:
            set: Vertices remaining after k-core filtering
        """
        # Track degree for each vertex
        vertex_degrees = {v: len(adjacency_map[v]) for v in vertex_list}
        active_vertices = set(vertex_list)
        removal_queue = deque()

        # Initialize queue with low-degree vertices
        for vertex in vertex_list:
            if vertex_degrees[vertex] < min_degree:
                removal_queue.append(vertex)

        # Iteratively remove vertices and update neighbors
        while removal_queue:
            current = removal_queue.popleft()

            if current not in active_vertices:
                continue

            active_vertices.remove(current)

            # Update degrees of neighboring vertices
            for neighbor in adjacency_map[current]:
                if neighbor in active_vertices:
                    vertex_degrees[neighbor] -= 1
                    if vertex_degrees[neighbor] < min_degree:
                        removal_queue.append(neighbor)

        return active_vertices

    @staticmethod
    def _count_components(vertex_set, connection_map):
        """
        Count connected components using breadth-first search

        Args:
            vertex_set: Set of vertices to analyze
            connection_map: Adjacency information

        Returns:
            int: Number of connected components
        """
        if not vertex_set:
            return 0

        visited_vertices = set()
        component_count = 0

        for start_vertex in vertex_set:
            if start_vertex in visited_vertices:
                continue

            # Start new component exploration
            component_count += 1
            exploration_queue = deque([start_vertex])
            visited_vertices.add(start_vertex)

            # BFS traversal for current component
            while exploration_queue:
                current_node = exploration_queue.popleft()

                for connected_node in connection_map[current_node]:
                    if (connected_node in vertex_set and
                        connected_node not in visited_vertices):
                        visited_vertices.add(connected_node)
                        exploration_queue.append(connected_node)

        return component_count

    ################################################################################
    # All auxiliary methods are contained within kCoreBaseStructuralDiversity class
    ################################################################################


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
