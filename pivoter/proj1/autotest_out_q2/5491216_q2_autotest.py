#!/usr/bin/env python3
# Auto-generated for 5491216

STUDENT_ID = "5491216"
STUDENT_NAME = "Haowei Peng"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import defaultdict
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
        diversity_measures = []
        graph_size = G.vertex_num
        
        for node_idx in range(graph_size):
            neighbor_list = G.adj_list[node_idx]
            structural_value = kCoreBaseStructuralDiversity._analyze_vertex_neighborhood(
                G, node_idx, neighbor_list, k)
            diversity_measures.append(structural_value)
        
        #print(diversity_measures)
        return diversity_measures

    @staticmethod
    def _analyze_vertex_neighborhood(graph_obj, target_node, adjacency_data, threshold_k):
        """Compute structural diversity for a specific vertex based on its neighborhood"""
        if not adjacency_data:
            return 0
            
        neighborhood_graph = kCoreBaseStructuralDiversity._construct_neighborhood_subgraph(
            graph_obj, adjacency_data)
        
        component_count = kCoreBaseStructuralDiversity._determine_k_core_components(
            neighborhood_graph, threshold_k)
        
        return component_count

    @staticmethod
    def _construct_neighborhood_subgraph(original_graph, vertex_neighbors):
        """Build induced subgraph from neighborhood vertices"""
        subgraph_edges = defaultdict(set)
        neighbor_lookup = set(vertex_neighbors)
        
        for current_node in vertex_neighbors:
            adjacent_nodes = original_graph.adj_list[current_node]
            for potential_neighbor in adjacent_nodes:
                if potential_neighbor in neighbor_lookup:
                    subgraph_edges[current_node].add(potential_neighbor)
                    
        return {node: list(edges) for node, edges in subgraph_edges.items()}

    @staticmethod 
    def _determine_k_core_components(edge_structure, min_connectivity):
        """Calculate number of k-core connected components"""
        if not edge_structure:
            return 0
            
        surviving_nodes = kCoreBaseStructuralDiversity._apply_k_core_filter(
            edge_structure, min_connectivity)
        
        if not surviving_nodes:
            return 0
            
        return kCoreBaseStructuralDiversity._enumerate_components(
            edge_structure, surviving_nodes)

    @staticmethod
    def _apply_k_core_filter(adjacency_structure, degree_bound):
        """Filter vertices maintaining k-core property through iterative removal"""
        node_connectivity = {vertex: len(connections) 
                           for vertex, connections in adjacency_structure.items()}
        
        active_vertex_set = set(adjacency_structure.keys())
        removal_candidates = []
        
        # Initialize removal queue with low-degree vertices
        for vertex in node_connectivity:
            if node_connectivity[vertex] < degree_bound:
                removal_candidates.append(vertex)
        
        # Iteratively remove vertices until stability
        while removal_candidates:
            current_vertex = removal_candidates.pop()
            
            if current_vertex not in active_vertex_set:
                continue
                
            active_vertex_set.remove(current_vertex)
            
            # Update connectivity of remaining neighbors
            for neighbor_vertex in adjacency_structure[current_vertex]:
                if neighbor_vertex in active_vertex_set:
                    node_connectivity[neighbor_vertex] -= 1
                    if node_connectivity[neighbor_vertex] < degree_bound:
                        removal_candidates.append(neighbor_vertex)
        
        return active_vertex_set

    @staticmethod
    def _enumerate_components(graph_topology, valid_vertex_set):
        """Count connected components using depth-first exploration"""
        unexplored_vertices = set(valid_vertex_set)
        total_components = 0
        
        while unexplored_vertices:
            total_components += 1
            exploration_root = unexplored_vertices.pop()
            vertex_stack = [exploration_root]
            
            # DFS traversal of current component
            while vertex_stack:
                current_vertex = vertex_stack.pop()
                
                for adjacent_vertex in graph_topology.get(current_vertex, []):
                    if adjacent_vertex in unexplored_vertices:
                        unexplored_vertices.remove(adjacent_vertex)
                        vertex_stack.append(adjacent_vertex)
        
        return total_components

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
