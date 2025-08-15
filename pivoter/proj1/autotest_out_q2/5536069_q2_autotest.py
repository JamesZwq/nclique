#!/usr/bin/env python3
# Auto-generated for 5536069

STUDENT_ID = "5536069"
STUDENT_NAME = "Jinzhe Wang"

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
        
        # Find k-core of the entire graph
        kcore_vertices = kCoreBaseStructuralDiversity._find_kcore(G, k)
        
        if not kcore_vertices:
            return [0] * n
        
        # Create mapping from original vertex id to kcore subgraph id
        kcore_list = sorted(kcore_vertices)
        vertex_to_kcore_id = {v: i for i, v in enumerate(kcore_list)}
        
        # Build adjacency list for k-core subgraph
        kcore_adj = [[] for _ in range(len(kcore_list))]
        for u in kcore_list:
            u_id = vertex_to_kcore_id[u]
            for v in G.adj_list[u]:
                if v in kcore_vertices:
                    v_id = vertex_to_kcore_id[v]
                    kcore_adj[u_id].append(v_id)
        
        # For each vertex, compute structural diversity
        sd = [0] * n
        
        for orig_v in range(n):
            if orig_v not in kcore_vertices:
                sd[orig_v] = 0
                continue
            
            v_kcore_id = vertex_to_kcore_id[orig_v]
            
            # Get neighbors of v in the k-core
            neighbors = set(kcore_adj[v_kcore_id])
            
            if len(neighbors) == 0:
                sd[orig_v] = 0
                continue
            
            # Create neighbor-induced subgraph
            neighbor_list = sorted(neighbors)
            neighbor_to_id = {v: i for i, v in enumerate(neighbor_list)}
            
            neighbor_adj = [[] for _ in range(len(neighbor_list))]
            for u_kcore in neighbor_list:
                u_id = neighbor_to_id[u_kcore]
                for v_kcore in kcore_adj[u_kcore]:
                    if v_kcore in neighbors:
                        v_id = neighbor_to_id[v_kcore]
                        neighbor_adj[u_id].append(v_id)
            
            # Create a temporary graph object for the neighbor-induced subgraph
            class TempGraph:
                def __init__(self, adj_list):
                    self.vertex_num = len(adj_list)
                    self.adj_list = adj_list
            
            neighbor_graph = TempGraph(neighbor_adj)
            
            # Find all k-cores in the neighbor-induced subgraph
            neighbor_kcores = kCoreBaseStructuralDiversity._find_all_kcores(neighbor_graph, k)
            
            # Count the number of different k-cores
            sd[orig_v] = len(neighbor_kcores)
        
        return sd
    
    @staticmethod
    def _find_kcore(G, k):
        """
        Find k-core of graph G using iterative vertex removal.
        Returns set of vertices in the k-core.
        """
        # Initialize degrees
        degrees = [len(G.adj_list[v]) for v in range(G.vertex_num)]
        active = set(range(G.vertex_num))
        queue = deque()
        
        # Initialize queue with vertices having degree < k
        for v in range(G.vertex_num):
            if degrees[v] < k:
                queue.append(v)
        
        # Iteratively remove vertices with degree < k
        while queue:
            v = queue.popleft()
            if v not in active:
                continue
                
            active.remove(v)
            
            # Update neighbors' degrees
            for u in G.adj_list[v]:
                if u in active:
                    degrees[u] -= 1
                    if degrees[u] < k:
                        queue.append(u)
        
        return active
    
    @staticmethod
    def _find_all_kcores(G, k):
        """
        Find all maximal k-cores in graph G.
        Returns list of sets, each representing a connected k-core component.
        """
        kcore_vertices = kCoreBaseStructuralDiversity._find_kcore(G, k)
        
        if not kcore_vertices:
            return []
        
        # Find connected components in the k-core subgraph
        visited = set()
        components = []
        
        def dfs(start):
            component = set()
            stack = [start]
            
            while stack:
                v = stack.pop()
                if v in visited:
                    continue
                    
                visited.add(v)
                component.add(v)
                
                for u in G.adj_list[v]:
                    if u in kcore_vertices and u not in visited:
                        stack.append(u)
            
            return component
        
        for v in kcore_vertices:
            if v not in visited:
                component = dfs(v)
                if component:
                    components.append(component)
        
        return components


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
