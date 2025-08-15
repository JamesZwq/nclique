#!/usr/bin/env python3
# Auto-generated for 5496331

STUDENT_ID = "5496331"
STUDENT_NAME = "Zhaohan Wang"

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
        vertex_count = G.vertex_num
        diversity_scores = [0] * vertex_count
        
        # Compute structural diversity for each vertex
        for vertex_id in range(vertex_count):
            adjacent_vertices = G.adj_list[vertex_id]
            
            # Handle isolated vertices
            if not adjacent_vertices:
                diversity_scores[vertex_id] = 0
                continue
            
            # Construct induced subgraph from neighbors
            neighbor_lookup = set(adjacent_vertices)
            induced_graph = {}
            
            # Build adjacency structure for induced subgraph
            for node in adjacent_vertices:
                induced_graph[node] = []
            
            # Populate edges within neighbor set
            for node in adjacent_vertices:
                for connected_node in G.adj_list[node]:
                    if connected_node in neighbor_lookup and connected_node != node:
                        induced_graph[node].append(connected_node)
            
            # Calculate k-core components in induced subgraph
            diversity_scores[vertex_id] = kCoreBaseStructuralDiversity._compute_core_components(
                induced_graph, adjacent_vertices, k
            )
        
        return diversity_scores
    
    @staticmethod
    def _compute_core_components(graph_structure, node_list, target_k):
        """
        Determine number of k-core components in given graph structure
        """
        if not node_list:
            return 0
        
        # Special case: k=0 means count all connected components
        if target_k == 0:
            return kCoreBaseStructuralDiversity._find_component_count(graph_structure, node_list)
        
        # Execute k-core decomposition algorithm
        active_nodes = set(node_list)
        decomposition_ongoing = True
        
        while decomposition_ongoing:
            decomposition_ongoing = False
            nodes_to_eliminate = []
            
            for current_node in active_nodes:
                # Calculate effective degree within active subgraph
                effective_degree = len([n for n in graph_structure[current_node] if n in active_nodes])
                if effective_degree < target_k:
                    nodes_to_eliminate.append(current_node)
            
            if nodes_to_eliminate:
                decomposition_ongoing = True
                active_nodes.difference_update(nodes_to_eliminate)
        
        # Return component count in remaining k-core structure
        if not active_nodes:
            return 0
        
        return kCoreBaseStructuralDiversity._find_component_count(graph_structure, list(active_nodes))
    
    @staticmethod
    def _find_component_count(adjacency_map, vertex_set):
        """
        Use breadth-first traversal to count connected components
        """
        if not vertex_set:
            return 0
        
        explored_vertices = set()
        total_components = 0
        
        for root_vertex in vertex_set:
            if root_vertex not in explored_vertices:
                # Initialize new component exploration
                total_components += 1
                exploration_queue = deque([root_vertex])
                explored_vertices.add(root_vertex)
                
                while exploration_queue:
                    current_vertex = exploration_queue.popleft()
                    
                    for adjacent_vertex in adjacency_map[current_vertex]:
                        if adjacent_vertex not in explored_vertices and adjacent_vertex in vertex_set:
                            explored_vertices.add(adjacent_vertex)
                            exploration_queue.append(adjacent_vertex)
        
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
