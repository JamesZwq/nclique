#!/usr/bin/env python3
# Auto-generated for 5545694

STUDENT_ID = "5545694"
STUDENT_NAME = "Li Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass
    
    @staticmethod
    def process(G, k):
        """
        Compute k-core based structural diversity for all vertices in graph G.
        
        For each vertex v, compute τ_k(v) = number of k-cores in neighbor-induced subgraph G[N(v)]
        
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        n = G.vertex_num
        sd = []
        
        # Compute τ_k(v) for each vertex v
        for v in range(n):
            neighbors = G.adj_list[v]
            tau_k_v = kCoreBaseStructuralDiversity._compute_tau_k(neighbors, G, k)
            sd.append(tau_k_v)
        
        return sd
    
    @staticmethod
    def _compute_tau_k(neighbors, G, k):
        """
        Compute τ_k(v) for vertex v given its neighbors.
        
        Steps:
        1. Build neighbor-induced subgraph G[N(v)]
        2. Find k-cores by removing vertices with degree < k
        3. Count connected components (k-cores)
        
        Parameters
        ----------
        neighbors : list of int
            Neighbor vertex IDs of vertex v
        G : UndirectedUnweightedGraph
            Original graph
        k : int
            k-core parameter
            
        Returns
        -------
        int : τ_k(v) - number of k-cores in neighbor-induced subgraph
        """
        if len(neighbors) == 0:
            return 0
        
        # Build neighbor-induced subgraph
        neighbor_set = set(neighbors)
        neighbor_to_index = {neighbor: i for i, neighbor in enumerate(neighbors)}
        
        # Create adjacency list for neighbor-induced subgraph
        nbr_adj = [[] for _ in range(len(neighbors))]
        for i, u in enumerate(neighbors):
            for neighbor_of_u in G.adj_list[u]:
                if neighbor_of_u in neighbor_set:
                    j = neighbor_to_index[neighbor_of_u]
                    if i != j:  # Avoid self-loops
                        nbr_adj[i].append(j)
        
        # Count k-cores in neighbor-induced subgraph
        return kCoreBaseStructuralDiversity._count_k_cores(nbr_adj, k)
    
    @staticmethod
    def _count_k_cores(adj_list, k):
        """
        Count the number of k-cores in a graph.
        
        A k-core is a maximal subgraph where every vertex has degree >= k.
        Algorithm:
        1. Iteratively remove vertices with degree < k
        2. Count connected components in remaining graph
        
        Parameters
        ----------
        adj_list : list of list of int
            Adjacency list representation of the graph
        k : int
            k-core parameter
            
        Returns
        -------
        int : number of k-cores (connected components after pruning)
        """
        n = len(adj_list)
        if n == 0:
            return 0
        
        # Create mutable copy of adjacency lists
        graph = [list(neighbors) for neighbors in adj_list]
        active = [True] * n  # Track active vertices
        
        # Iteratively remove vertices with degree < k
        while True:
            to_remove = []
            for v in range(n):
                if active[v] and len(graph[v]) < k:
                    to_remove.append(v)
            
            if not to_remove:
                break
            
            # Remove vertices and update adjacency lists
            for v in to_remove:
                active[v] = False
                # Remove v from all its neighbors' adjacency lists
                for neighbor in list(graph[v]):  # Use list() to avoid modification during iteration
                    if active[neighbor] and v in graph[neighbor]:
                        graph[neighbor].remove(v)
                graph[v] = []
        
        # Count connected components in remaining graph
        remaining_vertices = [v for v in range(n) if active[v]]
        if not remaining_vertices:
            return 0
        
        visited = [False] * n
        components = 0
        
        for v in remaining_vertices:
            if not visited[v]:
                # BFS to find connected component
                queue = deque([v])
                visited[v] = True
                
                while queue:
                    current = queue.popleft()
                    for neighbor in graph[current]:
                        if not visited[neighbor] and active[neighbor]:
                            visited[neighbor] = True
                            queue.append(neighbor)
                
                components += 1
        
        return components

################################################################################
# You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
################################################################################


################################################################################
# Do not edit this code cell.
from itertools import combinations, accumulate
################################################################################

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
