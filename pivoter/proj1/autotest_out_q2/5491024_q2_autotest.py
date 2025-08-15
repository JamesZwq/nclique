#!/usr/bin/env python3
# Auto-generated for 5491024

STUDENT_ID = "5491024"
STUDENT_NAME = "Jiajun Tan"

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
      
        n = G.vertex_num
        sd = [0] * n
        
        # Compute k-core based structural diversity for each vertex
        for v in range(n):
            # Get neighbor set of vertex v
            neighbors = set(G.adj_list[v])
            
            # Handle edge case: if not enough neighbors for k-core
            if len(neighbors) < k:
                sd[v] = 0
                continue
            
            # Build neighbor-induced subgraph
            # Map original vertex IDs to subgraph IDs
            neighbor_list = list(neighbors)
            vertex_map = {orig: idx for idx, orig in enumerate(neighbor_list)}
            subgraph_size = len(neighbor_list)
            
            # Build adjacency list for neighbor-induced subgraph
            subgraph_adj = [[] for _ in range(subgraph_size)]
            for i, u in enumerate(neighbor_list):
                for neighbor in G.adj_list[u]:
                    if neighbor in vertex_map:
                        subgraph_adj[i].append(vertex_map[neighbor])
            
            # Find k-cores in the neighbor-induced subgraph
            k_core_vertices = kCoreBaseStructuralDiversity._find_k_core_vertices(subgraph_adj, k)
            
            # Count connected components in k-core vertices
            sd[v] = kCoreBaseStructuralDiversity._count_connected_components(subgraph_adj, k_core_vertices)
        
        return sd
    
    @staticmethod
    def _find_k_core_vertices(adj_list, k):
        """
        Find vertices that belong to k-core using iterative removal.
        Returns set of vertex indices that form k-core.
        """
        n = len(adj_list)
        # Track current degree of each vertex
        degrees = [len(adj_list[i]) for i in range(n)]
        # Track which vertices are still active
        active = [True] * n
        # Queue for vertices to be removed
        removal_queue = deque()
        
        # Initialize queue with vertices having degree < k
        for i in range(n):
            if degrees[i] < k:
                removal_queue.append(i)
        
        # Iteratively remove vertices with degree < k
        while removal_queue:
            u = removal_queue.popleft()
            if not active[u]:
                continue
            
            # Mark vertex as removed
            active[u] = False
            
            # Update degrees of neighbors
            for neighbor in adj_list[u]:
                if active[neighbor]:
                    degrees[neighbor] -= 1
                    # If neighbor's degree drops below k, add to removal queue
                    if degrees[neighbor] < k:
                        removal_queue.append(neighbor)
        
        # Return set of vertices that remain (form k-core)
        return {i for i in range(n) if active[i]}
    
    @staticmethod
    def _count_connected_components(adj_list, vertex_set):
        """
        Count connected components among given vertex set.
        Uses DFS to find connected components.
        """
        if not vertex_set:
            return 0
        
        visited = set()
        component_count = 0
        
        # Perform DFS for each unvisited vertex
        for start_vertex in vertex_set:
            if start_vertex in visited:
                continue
            
            # Start new connected component
            component_count += 1
            dfs_stack = [start_vertex]
            
            # DFS traversal
            while dfs_stack:
                current = dfs_stack.pop()
                if current in visited:
                    continue
                
                visited.add(current)
                
                # Add unvisited neighbors to stack
                for neighbor in adj_list[current]:
                    if neighbor in vertex_set and neighbor not in visited:
                        dfs_stack.append(neighbor)
        
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
