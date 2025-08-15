#!/usr/bin/env python3
# Auto-generated for 5567959

STUDENT_ID = "5567959"
STUDENT_NAME = "Ember Xue"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity:

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

        # Basic input validation
        n = G.vertex_num
        sd = [0] * n

        adj = G.adj_list
        if k < 0:
            raise ValueError("k must be non-negative")
        # Handle edge cases
        if n == 0:
            return []
        if n == 1:
            return [0]



        #  Calculate global coreness for all vertices
        global_coreness = kCoreBaseStructuralDiversity.calculate_kcore_peeling(
            adj, target_k=None
        )

        # Process each vertex to find its structural diversity
        for v in range(n):
            # Get neighbors of current vertex
            neighbor_list = adj[v]

            # Skip if no neighbors
            if not neighbor_list:
                continue

            # Filter neighbors that have coreness >= k
            valid_neighbors = [
                neighbor for neighbor in neighbor_list
                if neighbor < n and global_coreness.get(neighbor, 0) >= k
            ]

            # Skip if no valid neighbors
            if not valid_neighbors:
                continue

            # Build induced subgraph from neighbors
            neighbor_subgraph = kCoreBaseStructuralDiversity.build_induced_subgraph(
                neighbor_list, adj
            )

            #Apply k-core peeling to find remaining nodes
            remaining_nodes_after_peeling = kCoreBaseStructuralDiversity.calculate_kcore_peeling(
                neighbor_subgraph, target_k=k
            )

            # If no nodes remain after peeling, diversity is 0
            if not remaining_nodes_after_peeling:
                continue

            # Count connected components in remaining k-core nodes
            component_count = kCoreBaseStructuralDiversity.count_connected_components(
                neighbor_subgraph, remaining_nodes_after_peeling
            )
            sd[v] = component_count

        return sd

    @staticmethod
    def calculate_kcore_peeling(graph_adj, target_k=None):
        """
        Universal k-core peeling algorithm.
        """
        # Handle empty input
        if not graph_adj:
            return {} if target_k is None else []

        # Determine node set and neighbor access method
        if isinstance(graph_adj, dict):
            node_list = list(graph_adj.keys())
            get_neighbors = lambda vertex: graph_adj.get(vertex, [])
        else:
            node_list = list(range(len(graph_adj)))
            get_neighbors = lambda vertex: graph_adj[vertex] if vertex < len(graph_adj) else []

        # Initialize degree for each node
        deg = {vertex: len(get_neighbors(vertex)) for vertex in node_list}

        # Initialize processing queue and tracking sets
        processing_queue = deque()
        removed_nodes = set()
        coreness_values = {}

        # Add initial nodes to queue based on mode
        if target_k is not None:
            # Local mode: queue nodes with degree < target_k
            for vertex in node_list:
                if deg[vertex] < target_k:
                    processing_queue.append(vertex)
        else:
            # Global mode: queue all nodes for coreness calculation
            for vertex in node_list:
                processing_queue.append(vertex)

        # Main peeling loop
        while processing_queue:
            curr_v = processing_queue.popleft()

            # Skip if already processed
            if curr_v in removed_nodes:
                continue

            removed_nodes.add(curr_v)

            # Record coreness for global mode
            if target_k is None:
                coreness_values[curr_v] = deg[curr_v]

            # Update degrees of neighbors after removing current vertex
            for neighbor in get_neighbors(curr_v):
                if neighbor in deg and neighbor not in removed_nodes:
                    deg[neighbor] -= 1

                    # Add neighbor to queue if conditions are met
                    if target_k is not None:
                        # Local mode: add if degree drops below threshold
                        if deg[neighbor] < target_k:
                            processing_queue.append(neighbor)
                    else:
                        # Global mode: add if degree drops below current vertex's degree
                        if deg[neighbor] < coreness_values[curr_v]:
                            processing_queue.append(neighbor)

        # Return results based on mode
        if target_k is None:
            # Assign final coreness to unprocessed nodes
            for vertex in node_list:
                if vertex not in coreness_values:
                    coreness_values[vertex] = deg.get(vertex, 0)
            return coreness_values
        else:
            # Return nodes that survived k-core peeling
            surviving_nodes = [vertex for vertex in node_list if vertex not in removed_nodes]
            return surviving_nodes


    @staticmethod
    def build_induced_subgraph(node_subset, adj):
        """
        Build an induced subgraph containing only edges within the node subset.
        Optimized with a boolean marker for O(1) membership checks.
        """
        if not node_subset:
            return {}

        n = len(adj)  # total number of nodes in the graph

        # Mark all nodes in the subset using a boolean array
        # This allows O(1) membership checking without hash lookups
        is_in_subset = [False] * n
        for u in node_subset:
            is_in_subset[u] = True

        subgraph_adj = {}

        # For each node in the subset, filter only neighbors that also belong to the subset
        for u in node_subset:
            if u >= n:
                continue  # safety check for out-of-bound nodes
            subgraph_adj[u] = [
                nb for nb in adj[u] if nb < n and is_in_subset[nb]
            ]

        #Clear the marker to avoid polluting subsequent logic
        for u in node_subset:
            is_in_subset[u] = False

        return subgraph_adj



    @staticmethod
    def count_connected_components(subgraph_adj, target_nodes):
        """
        Count the number of connected components among target nodes.
        """
        if not target_nodes:
            return 0

        target_node_set = set(target_nodes)
        visited = set()
        component_count = 0

        # Process each unvisited node as potential component start
        for start_node in target_nodes:
            if start_node not in visited:
                component_count += 1

                # BFS to find all nodes in this component
                bfs_queue = deque([start_node])
                visited.add(start_node)

                while bfs_queue:
                    current_node = bfs_queue.popleft()

                    # Check all neighbors of current node
                    neighbors = subgraph_adj.get(current_node, [])
                    for neighbor in neighbors:
                        if neighbor in target_node_set and neighbor not in visited:
                            visited.add(neighbor)
                            bfs_queue.append(neighbor)

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
