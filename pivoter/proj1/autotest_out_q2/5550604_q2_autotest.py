#!/usr/bin/env python3
# Auto-generated for 5550604

STUDENT_ID = "5550604"
STUDENT_NAME = "Pengpai Zhang"

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
        
        for v in range(n):
            # Get neighbor-induced subgraph for vertex v
            neighbors = set(G.adj_list[v])
            if len(neighbors) == 0:
                sd[v] = 0
                continue
                
            # Create mapping from original vertex id to local id in subgraph
            neighbor_list = list(neighbors)
            local_to_global = {i: neighbor_list[i] for i in range(len(neighbor_list))}
            global_to_local = {neighbor_list[i]: i for i in range(len(neighbor_list))}
            
            # Build adjacency list for neighbor-induced subgraph
            sub_adj = [[] for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list):
                for w in G.adj_list[u]:
                    if w in neighbors:  # w is also a neighbor of v
                        sub_adj[i].append(global_to_local[w])
            
            # Find k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity._find_k_cores(sub_adj, k)
            sd[v] = len(k_cores)
        
        return sd
    
    @staticmethod
    def _find_k_cores(adj_list, k):
        n = len(adj_list)
        if n == 0:
            return []
        
        # Find the maximal k-core using core decomposition
        degrees = [len(adj_list[i]) for i in range(n)]
        removed = [False] * n
        queue = deque()
        
        # Initialize queue with vertices having degree < k
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)
                removed[i] = True
        
        # Remove vertices with degree < k iteratively
        while queue:
            u = queue.popleft()
            for v in adj_list[u]:
                if not removed[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)
                        removed[v] = True
        
        # Collect remaining vertices (vertices in k-core)
        k_core_vertices = []
        for i in range(n):
            if not removed[i]:
                k_core_vertices.append(i)
        
        if len(k_core_vertices) == 0:
            return []
        
        # Build subgraph of k-core vertices
        vertex_to_idx = {v: i for i, v in enumerate(k_core_vertices)}
        k_core_adj = [[] for _ in range(len(k_core_vertices))]
        
        for i, v in enumerate(k_core_vertices):
            for u in adj_list[v]:
                if not removed[u]:
                    k_core_adj[i].append(vertex_to_idx[u])
        
        # Find all connected components in the k-core subgraph
        visited = [False] * len(k_core_vertices)
        components = []
        
        def dfs(start, component):
            stack = [start]
            while stack:
                node = stack.pop()
                if visited[node]:
                    continue
                visited[node] = True
                component.append(k_core_vertices[node])  # Add original vertex id
                for neighbor in k_core_adj[node]:
                    if not visited[neighbor]:
                        stack.append(neighbor)
        
        for i in range(len(k_core_vertices)):
            if not visited[i]:
                component = []
                dfs(i, component)
                if component:
                    components.append(set(component))
        
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
