#!/usr/bin/env python3
# Auto-generated for 5524512

STUDENT_ID = "5524512"
STUDENT_NAME = "Wenbo Li"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        
        # For each vertex, compute its k-core based structural diversity
        for v in range(n):
            neighbors = list(G.adj_list[v])
            if not neighbors:
                continue
            
            # Build neighbor-induced subgraph
            # Create mapping from original node ID to local index
            neighbor_map = {}
            for i, node in enumerate(neighbors):
                neighbor_map[node] = i
            
            # Initialize adjacency list for subgraph
            subgraph_adj = []
            for i in range(len(neighbors)):
                subgraph_adj.append([])
            
            # Build edges between neighbors
            for i, u in enumerate(neighbors):
                for neighbor in G.adj_list[u]:
                    if neighbor in neighbor_map:
                        local_id = neighbor_map[neighbor]
                        subgraph_adj[i].append(local_id)
            
            # Perform k-core decomposition on neighbor-induced subgraph
            remaining = kCoreBaseStructuralDiversity.decompose_k_core(subgraph_adj, k)
            
            # Count connected components in remaining vertices
            sd[v] = kCoreBaseStructuralDiversity.count_connected_components(subgraph_adj, remaining)
        
        return sd

    @staticmethod
    def decompose_k_core(adj_list, k):
        """Extract k-core vertices from graph using iterative removal"""
        n = len(adj_list)
        if n == 0:
            return set()
        
        # Calculate initial degree for each vertex
        degrees = []
        for i in range(n):
            degree = len(adj_list[i])
            degrees.append(degree)
        
        # Track which vertices are removed
        removed = []
        for i in range(n):
            removed.append(False)
        
        # Initialize queue with low-degree vertices
        queue = deque()
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)
                removed[i] = True
        
        # Iteratively remove vertices
        while queue:
            vertex = queue.popleft()
            for neighbor in adj_list[vertex]:
                if not removed[neighbor]:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        queue.append(neighbor)
                        removed[neighbor] = True
        
        # Collect remaining vertices
        remaining_vertices = set()
        for i in range(n):
            if not removed[i]:
                remaining_vertices.add(i)
        
        return remaining_vertices

    @staticmethod
    def count_connected_components(adj_list, vertex_set):
        """Count connected components using BFS traversal"""
        if not vertex_set:
            return 0
        
        visited = set()
        component_count = 0
        
        # Try each vertex as a potential starting point
        for start_vertex in vertex_set:
            # Skip if already visited
            if start_vertex in visited:
                continue
                
            # Start BFS exploration for new component
            queue = deque()
            queue.append(start_vertex)
            visited.add(start_vertex)
            component_count += 1
            
            # Explore all vertices in this component
            while queue:
                current = queue.popleft()
                
                # Check all neighbors of current vertex
                for neighbor in adj_list[current]:
                    # Add neighbor if it's in vertex_set and not visited
                    if neighbor in vertex_set:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
        
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
