#!/usr/bin/env python3
# Auto-generated for 5503669

STUDENT_ID = "5503669"
STUDENT_NAME = "(Ryan) Borui Meng"

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
        tau = [0] * n
        
        core_numbers = kCoreBaseStructuralDiversity._compute_core_numbers(G)
        
        for v in range(n):
            tau[v] = kCoreBaseStructuralDiversity._compute_structural_diversity_optimized(G, v, k, core_numbers)
        
        return tau
    
    @staticmethod
    def _compute_core_numbers(G):
        """Compute core numbers using bucket-sort peeling algorithm"""
        n = G.vertex_num
        if n == 0:
            return []
        
        degrees = [len(G.adj_list[v]) for v in range(n)]
        core_numbers = [0] * n
        removed = [False] * n
        
        max_degree = max(degrees) if degrees else 0
        buckets = [[] for _ in range(max_degree + 1)]
        
        for v in range(n):
            buckets[degrees[v]].append(v)
        
        for current_k in range(max_degree + 1):
            while buckets[current_k]:
                v = buckets[current_k].pop()
                if removed[v]:
                    continue
                
                core_numbers[v] = current_k
                removed[v] = True
                
                for u in G.adj_list[v]:
                    if not removed[u]:
                        degrees[u] -= 1
                        if degrees[u] >= 0:
                            buckets[degrees[u]].append(u)
        
        return core_numbers
    
    @staticmethod
    def _compute_structural_diversity_optimized(G, v, k, core_numbers):
        """Compute structural diversity with pruning strategies"""
        neighbors = G.adj_list[v]
        
        if len(neighbors) < k:
            return 0
        
        valid_neighbors = [u for u in neighbors if core_numbers[u] >= k]
        if len(valid_neighbors) < k:
            return 0
        
        neighbor_adj = kCoreBaseStructuralDiversity._build_neighbor_subgraph_optimized(G, valid_neighbors)
        
        if not kCoreBaseStructuralDiversity._has_potential_k_core(neighbor_adj, k):
            return 0
        
        k_core_vertices = kCoreBaseStructuralDiversity._compute_k_core_optimized(neighbor_adj, k)
        
        return kCoreBaseStructuralDiversity._count_connected_components_optimized(k_core_vertices, neighbor_adj)
    
    @staticmethod
    def _build_neighbor_subgraph_optimized(G, valid_neighbors):
        """Build neighbor subgraph with filtered vertices"""
        neighbor_set = set(valid_neighbors)
        neighbor_adj = {}
        
        for u in valid_neighbors:
            adj_list = []
            for w in G.adj_list[u]:
                if w in neighbor_set:
                    adj_list.append(w)
            neighbor_adj[u] = adj_list
        
        return neighbor_adj
    
    @staticmethod
    def _has_potential_k_core(neighbor_adj, k):
        """Check if subgraph can contain k-core"""
        valid_count = 0
        for u in neighbor_adj:
            if len(neighbor_adj[u]) >= k:
                valid_count += 1
                if valid_count >= k:
                    return True
        return False
    
    @staticmethod
    def _compute_k_core_optimized(neighbor_adj, k):
        """k-core computation with early termination"""
        if not neighbor_adj:
            return set()
        
        active_vertices = set(neighbor_adj.keys())
        degrees = {u: len(neighbor_adj[u]) for u in neighbor_adj}
        
        if all(deg < k for deg in degrees.values()):
            return set()
        
        to_remove = deque()
        for u in active_vertices:
            if degrees[u] < k:
                to_remove.append(u)
        
        while to_remove:
            u = to_remove.popleft()
            
            if u not in active_vertices:
                continue
            
            active_vertices.remove(u)
            
            if len(active_vertices) < k:
                return set()
            
            for w in neighbor_adj[u]:
                if w in active_vertices:
                    degrees[w] -= 1
                    if degrees[w] < k:
                        to_remove.append(w)
        
        return active_vertices
    
    @staticmethod
    def _count_connected_components_optimized(vertices, neighbor_adj):
        """Count connected components using DFS"""
        if not vertices:
            return 0
        
        visited = set()
        component_count = 0
        
        for start_vertex in vertices:
            if start_vertex not in visited:
                component_count += 1
                
                stack = [start_vertex]
                while stack:
                    current = stack.pop()
                    
                    if current not in visited:
                        visited.add(current)
                        
                        for neighbor in neighbor_adj[current]:
                            if neighbor in vertices and neighbor not in visited:
                                stack.append(neighbor)
        
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
