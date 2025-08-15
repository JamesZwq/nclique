#!/usr/bin/env python3
# Auto-generated for 5598414

STUDENT_ID = "5598414"
STUDENT_NAME = "Yuzhou Zhu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
from itertools import combinations
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core based structural diversity for all vertices in G.
        
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
        
        # Precompute core values for all vertices
        core_values = kCoreBaseStructuralDiversity.core_decomposition(G)
        
        for v in range(n):
            # Get neighbor-induced subgraph vertices
            neighbors = G.adj_list[v]
            if not neighbors:
                continue
                
            # Compute k-core of neighbor-induced subgraph
            k_core = kCoreBaseStructuralDiversity.compute_k_core(G, v, neighbors, k, core_values)
            
            # Count connected components in k-core
            tau[v] = kCoreBaseStructuralDiversity.count_components(G, k_core)
            
        return tau
    
    @staticmethod
    def core_decomposition(G):
        """
        Implementation of Batagelj and Zaversnik's O(m) algorithm.
        """
        n = G.vertex_num
        degrees = [len(G.adj_list[v]) for v in range(n)]
        max_degree = max(degrees) if degrees else 0
        
        # Initialize data structures
        core = [0] * n
        bin_counts = [0] * (max_degree + 2)
        pos = [0] * n
        vert = [0] * n
        
        # Count vertices of each degree
        for v in range(n):
            bin_counts[degrees[v]] += 1
            
        # Set starting positions of bins
        start = 0
        for d in range(max_degree + 1):
            num = bin_counts[d]
            bin_counts[d] = start
            start += num
            
        # Sort vertices by degree
        for v in range(n):
            pos[v] = bin_counts[degrees[v]]
            vert[pos[v]] = v
            bin_counts[degrees[v]] += 1
            
        # Reset bin_counts to correct positions
        for d in range(max_degree, 0, -1):
            bin_counts[d] = bin_counts[d-1]
        bin_counts[0] = 0
        
        # Core decomposition
        for i in range(n):
            v = vert[i]
            core[v] = degrees[v]
            
            for u in G.adj_list[v]:
                if degrees[u] > degrees[v]:
                    du = degrees[u]
                    pu = pos[u]
                    pw = bin_counts[du]
                    w = vert[pw]
                    
                    if u != w:
                        pos[u] = pw
                        pos[w] = pu
                        vert[pu] = w
                        vert[pw] = u
                        
                    bin_counts[du] += 1
                    degrees[u] -= 1
                    
        return core
    
    @staticmethod
    def compute_k_core(G, center_v, vertices, k, core_values):
        """
        Implementation of k-core computation for neighbor-induced subgraph.
        """
        if not vertices:
            return set()
            
        # Create subgraph with only vertices that could be in k-core
        subgraph = defaultdict(set)
        degrees = {}
        remaining = set()
        
        for v in vertices:
            if core_values[v] <= k:
                continue  # Can't be in k-core
            remaining.add(v)
            
        # Build subgraph and initial degrees
        for v in remaining:
            neighbors_in_subgraph = set()
            for u in G.adj_list[v]:
                if u in remaining:
                    neighbors_in_subgraph.add(u)
            subgraph[v] = neighbors_in_subgraph
            degrees[v] = len(neighbors_in_subgraph)
        
        # Initialize queue with vertices having degree < k
        queue = deque()
        for v in remaining:
            if degrees[v] < k:
                queue.append(v)
                
        # Process queue - iterative removal
        while queue:
            v = queue.popleft()
            if v not in remaining:
                continue
                
            # Remove v from subgraph
            remaining.remove(v)
            
            # Update neighbors' degrees
            for u in subgraph[v]:
                if u in remaining:
                    degrees[u] -= 1
                    if degrees[u] < k:
                        queue.append(u)
        
        return remaining
    
    @staticmethod
    def count_components(G, vertices):
        """
        Count connected components in the given set of vertices.
        Uses BFS to traverse components.
        """
        if not vertices:
            return 0
            
        visited = set()
        component_count = 0
        
        for v in vertices:
            if v not in visited:
                # Start new component
                component_count += 1
                queue = deque([v])
                visited.add(v)
                
                while queue:
                    current = queue.popleft()
                    for neighbor in G.adj_list[current]:
                        if neighbor in vertices and neighbor not in visited:
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
