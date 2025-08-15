#!/usr/bin/env python3
# Auto-generated for 5485499

STUDENT_ID = "5485499"
STUDENT_NAME = "(Carol) Xue Zhang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################
class kCoreBaseStructuralDiversity:
    
    @staticmethod
    def process(graph, k_threshold):
        diversity_counts = []
        total_vertices = graph.vertex_num
        
        for vertex in range(total_vertices):
            neighbor_subgraph = kCoreBaseStructuralDiversity._build_neighbor_subgraph(graph, vertex)
            core_components = kCoreBaseStructuralDiversity._analyze_core_components(
                neighbor_subgraph, 
                k_threshold
            )
            diversity_counts.append(core_components)
            
        return diversity_counts

    @staticmethod
    def _build_neighbor_subgraph(graph, center_vertex):
        """Constructs the subgraph induced by a vertex's neighbors."""
        neighbors = graph.adj_list[center_vertex]
        neighbor_set = set(neighbors)
        subgraph = defaultdict(list)
        
        for node in neighbors:
            connected_neighbors = []
            for neighbor in graph.adj_list[node]:
                if neighbor in neighbor_set:
                    connected_neighbors.append(neighbor)
            if connected_neighbors:
                subgraph[node] = connected_neighbors
                    
        return subgraph

    @staticmethod
    def _analyze_core_components(subgraph, min_core_size):
        """Identifies and counts k-core components in a subgraph."""
        if not subgraph:
            return 0

        degree_map = {}
        for node in subgraph:
            degree_map[node] = len(subgraph[node])

        core_nodes = kCoreBaseStructuralDiversity._perform_core_decomposition(
            subgraph, 
            degree_map, 
            min_core_size
        )
        
        return kCoreBaseStructuralDiversity._count_connected_components(
            subgraph, 
            core_nodes
        )

    @staticmethod
    def _perform_core_decomposition(subgraph, degree_map, threshold):
        """Identifies nodes belonging to the k-core."""
        removal_queue = deque()
        pruned_nodes = set()

        for node in degree_map:
            if degree_map[node] < threshold:
                removal_queue.append(node)
                pruned_nodes.add(node)

        while removal_queue:
            current_node = removal_queue.popleft()
            
            for neighbor in subgraph[current_node]:
                if neighbor not in pruned_nodes:
                    degree_map[neighbor] -= 1
                    if degree_map[neighbor] < threshold:
                        pruned_nodes.add(neighbor)
                        removal_queue.append(neighbor)

        core_nodes = set()
        for node in subgraph:
            if node not in pruned_nodes:
                core_nodes.add(node)

        return core_nodes

    @staticmethod
    def _count_connected_components(subgraph, component_nodes):
        """Counts connected components among core nodes."""
        visited_nodes = set()
        component_count = 0

        for node in component_nodes:
            if node not in visited_nodes:
                component_count += 1
                kCoreBaseStructuralDiversity._explore_component(
                    subgraph, 
                    node, 
                    visited_nodes, 
                    component_nodes
                )

        return component_count

    @staticmethod
    def _explore_component(subgraph, start_node, visited_nodes, valid_nodes):
        """Performs BFS to mark all nodes in a connected component."""
        node_queue = deque()
        node_queue.append(start_node)
        visited_nodes.add(start_node)
        
        while node_queue:
            current_node = node_queue.popleft()
            
            for neighbor in subgraph[current_node]:
                if neighbor in valid_nodes and neighbor not in visited_nodes:
                    visited_nodes.add(neighbor)
                    node_queue.append(neighbor)

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
