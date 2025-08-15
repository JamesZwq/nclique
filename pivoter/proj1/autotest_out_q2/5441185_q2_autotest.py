#!/usr/bin/env python3
# Auto-generated for 5441185

STUDENT_ID = "5441185"
STUDENT_NAME = "Peter Lin"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core based structural diversity for all vertices
        
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
        
        # For each vertex, compute its k-core based structural diversity
        for v in range(n):
            tau[v] = kCoreBaseStructuralDiversity._compute_tau_k(G, v, k)
        
        return tau
    
    @staticmethod
    def _compute_tau_k(G, v, k):
        """
        Compute τ_k(v) for a single vertex v
        
        Args:
            G: Graph
            v: vertex
            k: k-core parameter
        Returns:
            Number of k-cores in neighbor-induced subgraph of v
        """
        # Get neighbors of v
        neighbors = G.adj_list[v]
        neighbor_set = set(neighbors)
        
        if len(neighbors) == 0:
            return 0
        
        # Build neighbor-induced subgraph
        # Map neighbor vertices to 0-indexed positions
        neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbors)}
        neighbor_adj = [[] for _ in range(len(neighbors))]
        
        # Add edges between neighbors
        for i, u in enumerate(neighbors):
            for neighbor_of_u in G.adj_list[u]:
                if neighbor_of_u in neighbor_set and neighbor_of_u != u:
                    j = neighbor_to_idx[neighbor_of_u]
                    if j not in neighbor_adj[i]:  # Avoid duplicates
                        neighbor_adj[i].append(j)
        
        # Find k-cores in the neighbor-induced subgraph
        return kCoreBaseStructuralDiversity._count_k_cores(neighbor_adj, k)
    
    @staticmethod
    def _count_k_cores(adj_list, k):
        """
        Count the number of k-cores in a graph represented by adjacency list
        
        Args:
            adj_list: adjacency list representation
            k: k-core parameter
        Returns:
            Number of k-cores (connected components where each vertex has degree >= k)
        """
        n = len(adj_list)
        if n == 0:
            return 0
        
        # Create a copy of adjacency list for modification
        adj_copy = [list(neighbors) for neighbors in adj_list]
        degree = [len(neighbors) for neighbors in adj_copy]
        removed = [False] * n
        
        # Remove vertices with degree < k iteratively
        changed = True
        while changed:
            changed = False
            for u in range(n):
                if not removed[u] and degree[u] < k:
                    # Remove vertex u
                    removed[u] = True
                    changed = True
                    
                    # Update degrees of its neighbors
                    for v in adj_copy[u]:
                        if not removed[v]:
                            degree[v] -= 1
                            # Remove u from v's adjacency list
                            adj_copy[v].remove(u)
        
        # Count connected components in remaining graph
        visited = [False] * n
        num_components = 0
        
        for u in range(n):
            if not removed[u] and not visited[u]:
                # Start BFS/DFS from this vertex
                kCoreBaseStructuralDiversity._dfs(u, adj_copy, visited, removed)
                num_components += 1
        
        return num_components
    
    @staticmethod
    def _dfs(start, adj_list, visited, removed):
        """
        DFS to mark all vertices in the same connected component
        """
        stack = [start]
        visited[start] = True
        
        while stack:
            u = stack.pop()
            for v in adj_list[u]:
                if not removed[v] and not visited[v]:
                    visited[v] = True
                    stack.append(v)

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
