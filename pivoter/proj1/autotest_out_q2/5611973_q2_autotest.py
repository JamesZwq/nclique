#!/usr/bin/env python3
# Auto-generated for 5611973

STUDENT_ID = "5611973"
STUDENT_NAME = "Ruicheng Luo"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
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
        res = [0] * n
        
        for v in range(n):
            # Get neighbors of v
            neighbors = G.adj_list[v]
            
            if len(neighbors) == 0:
                res[v] = 0
                continue
            
            # Build neighbor-induced subgraph
            neighbor_set = set(neighbors)
            subgraph_adj = defaultdict(list)
            
            # Build adjacency list
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbor_set and w != u:
                        subgraph_adj[u].append(w)
            
            # Using another function to find k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity._find_k_cores(subgraph_adj, neighbors, k)
            res[v] = len(k_cores)
        
        return res
    
    @staticmethod
    def _find_k_cores(adj_list, vertices, k):
        """
        Find all k-cores in the given graph.
        Returns a list of k-cores, where each k-core is a set of vertices.
        """
        # Create a copy
        adj = defaultdict(set)
        for v in vertices:
            if v in adj_list:
                adj[v] = set(adj_list[v])
        
        # Compute degrees
        degree = {}
        for v in vertices:
            if v in adj:
                degree[v] = len(adj[v]) 
            else:
                degree[v] = 0
        
        # Remove vertices with degree < k iteratively
        queue = deque()
        removed = set()
        
        # Initial queue with vertices which have degree < k 
        for v in vertices:
            if degree[v] < k:
                queue.append(v)
                removed.add(v)
        
        # Iteratively remove vertices
        while queue:
            v = queue.popleft()
            
            # Update degrees
            if v in adj:
                for u in adj[v]:
                    if u not in removed:
                        degree[u] -= 1
                        if degree[u] < k:
                            queue.append(u)
                            removed.add(u)
        
        k_core_vertices = set(vertices) - removed
        
        if not k_core_vertices:
            return []
        
        # Find connected component in k-core subgraph
        k_cores = []
        visited = set()
        
        for start in k_core_vertices:
            if start not in visited:
                # Use BFS to find connected component
                nodes = set()
                queue = deque([start])
                visited.add(start)
                nodes.add(start)
                
                while queue:
                    v = queue.popleft()
                    if v in adj:
                        for u in adj[v]:
                            if u in k_core_vertices and u not in visited:
                                visited.add(u)
                                queue.append(u)
                                nodes.add(u)
                
                # Verify this component is indeed a k-core (all vertices should have degree >= k)
                is_k_core = True
                for v in nodes:
                    nodes_degree = 0
                    if v in adj:
                        for u in adj[v]:
                            if u in nodes:
                                nodes_degree += 1
                    if nodes_degree < k:
                        is_k_core = False
                        break
                
                if is_k_core:
                    k_cores.append(nodes)
        
        return k_cores


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
