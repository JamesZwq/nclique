#!/usr/bin/env python3
# Auto-generated for 5463005

STUDENT_ID = "5463005"
STUDENT_NAME = "Tianxu Wang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
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
        num_vertices = G.vertex_num
        structural_diversity_values = [0] * num_vertices
        
        for current_vertex_id in range(num_vertices):
            # Retrieve neighbors of the current vertex
            vertex_neighbors = G.adj_list[current_vertex_id]
            if not vertex_neighbors:  # If no neighbors, diversity is 0
                structural_diversity_values[current_vertex_id] = 0
                continue
                
            # Construct a subgraph induced by the neighbors
            neighbor_set = set(vertex_neighbors)
            
            # Build adjacency list for this neighbor-induced subgraph
            neighbor_subgraph_adj = {}
            for nbr_u in vertex_neighbors:
                neighbor_subgraph_adj[nbr_u] = []
                for nbr_w in G.adj_list[nbr_u]:
                    if nbr_w in neighbor_set:
                        neighbor_subgraph_adj[nbr_u].append(nbr_w)
            
            # Calculate k-core based structural diversity for the current vertex
            structural_diversity_values[current_vertex_id] = kCoreBaseStructuralDiversity._compute_k_core_components(neighbor_subgraph_adj, k)
        
        return structural_diversity_values

    @staticmethod
    def _compute_k_core_components(graph_adjacency_dict, k_value):
        """
        Calculates the number of connected components within the k-core of a given graph.
        This function now also handles the connected component counting directly.
        """
        if not graph_adjacency_dict or k_value < 0:
            return 0
            
        if k_value == 0:
            # Handle k=0 case: count connected components in the original graph_adjacency_dict
            visited_for_cc = set()
            components_count_k0 = 0
            for node_k0 in graph_adjacency_dict:
                if node_k0 not in visited_for_cc:
                    bfs_queue_k0 = deque([node_k0])
                    visited_for_cc.add(node_k0)
                    while bfs_queue_k0:
                        current_node_k0 = bfs_queue_k0.popleft()
                        for neighbor_k0 in graph_adjacency_dict[current_node_k0]:
                            if neighbor_k0 not in visited_for_cc:
                                visited_for_cc.add(neighbor_k0)
                                bfs_queue_k0.append(neighbor_k0)
                    components_count_k0 += 1
            return components_count_k0
            
        # Create a mutable copy of the graph for k-core decomposition
        current_graph_state = {node: set(neighbors_list) for node, neighbors_list in graph_adjacency_dict.items()}
        
        has_changed = True
        while has_changed:
            has_changed = False
            nodes_to_prune = []
            
            # Identify nodes with degree less than k_value
            for node_id in current_graph_state:
                if len(current_graph_state[node_id]) < k_value:
                    nodes_to_prune.append(node_id)
            
            # Remove identified nodes and their incident edges
            for node_to_remove in nodes_to_prune:
                if node_to_remove in current_graph_state: # Ensure node still exists before processing
                    for adjacent_node in current_graph_state[node_to_remove]:
                        if adjacent_node in current_graph_state: # Ensure adjacent node still exists
                            current_graph_state[adjacent_node].discard(node_to_remove)
                    del current_graph_state[node_to_remove]
                    has_changed = True
        
        if not current_graph_state: # If graph is empty after k-core decomposition
            return 0
            
        # Count connected components in the resulting k-core
        nodes_visited_k_core = set()
        component_count_k_core = 0
        
        for node_in_k_core in current_graph_state:
            if node_in_k_core not in nodes_visited_k_core:
                bfs_queue_k_core = deque([node_in_k_core])
                nodes_visited_k_core.add(node_in_k_core)
                
                while bfs_queue_k_core:
                    current_node_bfs_k_core = bfs_queue_k_core.popleft()
                    for neighbor_node_k_core in current_graph_state[current_node_bfs_k_core]:
                        if neighbor_node_k_core not in nodes_visited_k_core:
                            nodes_visited_k_core.add(neighbor_node_k_core)
                            bfs_queue_k_core.append(neighbor_node_k_core)
                
                component_count_k_core += 1
                
        return component_count_k_core

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
