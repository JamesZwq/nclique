#!/usr/bin/env python3
# Auto-generated for 5512113

STUDENT_ID = "5512113"
STUDENT_NAME = "Zengqi Hao"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _find_connected_clusters(adjacency, retained_nodes, min_connections):

        vertex_count = len(adjacency)
        processed = [False] * vertex_count
        cluster_count = 0

        # Find all connected components via BFS
        for node_idx in range(vertex_count):
            if retained_nodes[node_idx] and not processed[node_idx]:
                # New component found
                cluster_count += 1
                search_queue = deque([node_idx])
                processed[node_idx] = True

                # Process all nodes in this component
                while search_queue:
                    current = search_queue.popleft()
                    for linked_node in adjacency[current]:
                        if retained_nodes[linked_node] and not processed[linked_node]:
                            processed[linked_node] = True
                            search_queue.append(linked_node)

        return cluster_count

    @staticmethod
    def _count_k_cores(connection_map, threshold):
        """
        Count the number of k-cores in a graph represented by adjacency list
        """
        total_nodes = len(connection_map)
        if total_nodes == 0:
            return 0

        # Calculate initial node degrees
        connection_counts = []
        for node_links in connection_map:
            connection_counts.append(len(node_links))

        # Track nodes that stay in the k-core
        remaining_nodes = [True] * total_nodes
        removal_queue = deque()

        # Identify nodes with degree below threshold
        for node_idx in range(total_nodes):
            if connection_counts[node_idx] < threshold:
                removal_queue.append(node_idx)
                remaining_nodes[node_idx] = False

        # Iteratively remove low-degree nodes and update neighbors
        while removal_queue:
            removed = removal_queue.popleft()
            # Update neighbors of the removed node
            for adjacent in connection_map[removed]:
                if remaining_nodes[adjacent]:
                    connection_counts[adjacent] -= 1
                    if connection_counts[adjacent] < threshold:
                        removal_queue.append(adjacent)
                        remaining_nodes[adjacent] = False

        # Find connected components in the filtered graph
        return kCoreBaseStructuralDiversity._find_connected_clusters(connection_map, remaining_nodes, threshold)

    @staticmethod
    def process(graph, min_degree):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        total_vertices = graph.vertex_num
        diversity_values = [0] * total_vertices

        # Process each vertex individually
        for vertex_idx in range(total_vertices):
            # Extract neighborhood
            adjacent_vertices = set(graph.adj_list[vertex_idx])

            # Skip isolated vertices
            if not adjacent_vertices:
                diversity_values[vertex_idx] = 0
                continue

            # Convert to ordered list for consistent indexing
            neighbor_array = list(adjacent_vertices)

            # Create mapping between original and local indices
            id_mapping = {}
            for local_idx, original_id in enumerate(neighbor_array):
                id_mapping[original_id] = local_idx

            # Build the neighborhood subgraph
            subgraph_connections = [[] for _ in range(len(neighbor_array))]

            # Populate subgraph with valid edges
            for i, original_id in enumerate(neighbor_array):
                for connected_node in graph.adj_list[original_id]:
                    if connected_node in adjacent_vertices:
                        subgraph_connections[i].append(id_mapping[connected_node])

            # Calculate structural diversity score
            diversity_values[vertex_idx] = kCoreBaseStructuralDiversity._count_k_cores(subgraph_connections, min_degree)

        return diversity_values

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
