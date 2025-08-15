#!/usr/bin/env python3
# Auto-generated for 5485407

STUDENT_ID = "5485407"
STUDENT_NAME = "Xinglun Zhou"

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
        n = G.vertex_num
        sd = [0] * n
        
        # Process each vertex independently
        for vertex_id in range(n):
            # Extract neighbor information
            neighbors = G.adj_list[vertex_id]
            
            if len(neighbors) == 0:
                sd[vertex_id] = 0
                continue
            
            # Build neighbor subgraph representation
            neighbor_graph = kCoreBaseStructuralDiversity._build_neighbor_subgraph(
                G, neighbors
            )
            
            # Compute structural diversity
            sd[vertex_id] = kCoreBaseStructuralDiversity._count_k_core_components(
                neighbor_graph, k
            )
        
        return sd

    @staticmethod
    def _build_neighbor_subgraph(original_graph, neighbor_list):
        """
        Construct the induced subgraph from the neighbor set.
        Returns a dictionary representation of the subgraph.
        """
        # Create a set for quick lookup
        neighbor_set = set(neighbor_list)
        
        # Build adjacency representation for subgraph
        subgraph = {}
        
        for node in neighbor_list:
            connections = []
            # Check connections within the neighbor set
            for adjacent in original_graph.adj_list[node]:
                if adjacent in neighbor_set:
                    connections.append(adjacent)
            subgraph[node] = connections
        
        return subgraph
    
    @staticmethod
    def _count_k_core_components(subgraph, k):
        """
        Count the number of connected components in the k-core of the subgraph.
        """
        if not subgraph:
            return 0
        
        # Handle k=0 case directly
        if k == 0:
            return kCoreBaseStructuralDiversity._count_components(subgraph)
        
        # Extract k-core
        k_core = kCoreBaseStructuralDiversity._extract_k_core(subgraph, k)
        
        # Count components in k-core
        return kCoreBaseStructuralDiversity._count_components(k_core)
    
    @staticmethod
    def _extract_k_core(graph, k):
        """
        Extract the k-core from the given graph using iterative pruning.
        """
        # Create a working copy
        working_graph = {}
        degrees = {}
        
        for node, neighbors in graph.items():
            working_graph[node] = set(neighbors)
            degrees[node] = len(neighbors)
        
        # Iterative pruning process
        changed = True
        while changed:
            changed = False
            nodes_to_remove = []
            
            # Find nodes with degree < k
            for node, degree in degrees.items():
                if degree < k and node in working_graph:
                    nodes_to_remove.append(node)
            
            # Remove low-degree nodes
            for node in nodes_to_remove:
                if node in working_graph:
                    # Update neighbors' degrees
                    for neighbor in working_graph[node]:
                        if neighbor in degrees:
                            degrees[neighbor] -= 1
                    
                    # Remove the node
                    del working_graph[node]
                    del degrees[node]
                    
                    # Remove edges to this node
                    for other_node in working_graph:
                        if node in working_graph[other_node]:
                            working_graph[other_node].remove(node)
                    
                    changed = True
        
        # Convert back to list format
        result = {}
        for node, neighbors in working_graph.items():
            result[node] = list(neighbors)
        
        return result
    
    @staticmethod
    def _count_components(graph):
        """
        Count connected components using BFS traversal.
        """
        if not graph:
            return 0
        
        visited = set()
        component_count = 0
        
        for start_node in graph:
            if start_node not in visited:
                # Start BFS from unvisited node
                queue = deque([start_node])
                visited.add(start_node)
                
                while queue:
                    current = queue.popleft()
                    
                    # Explore neighbors
                    for neighbor in graph.get(current, []):
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
                
                component_count += 1
        
        return component_count

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
