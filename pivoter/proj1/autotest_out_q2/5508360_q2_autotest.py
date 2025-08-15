#!/usr/bin/env python3
# Auto-generated for 5508360

STUDENT_ID = "5508360"
STUDENT_NAME = "Dawei Hao"

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
        List[int] # τ_k(v) for all v
        """
        n = G.vertex_num
        sd = [0] * n
        
        # For each vertex, compute its k-core based structural diversity
        for v in range(n):
            # Get neighbors of v
            neighbors = set(G.adj_list[v])
            if len(neighbors) == 0:
                sd[v] = 0
                continue
            
            # Build neighbor-induced subgraph
            # Map from original vertex id to local id in subgraph
            vertex_map = {u: i for i, u in enumerate(sorted(neighbors))}
            reverse_map = {i: u for u, i in vertex_map.items()}
            subgraph_size = len(neighbors)
            
            # Build adjacency list for neighbor-induced subgraph
            subgraph_adj = [[] for _ in range(subgraph_size)]
            for u in neighbors:
                local_u = vertex_map[u]
                for w in G.adj_list[u]:
                    if w in neighbors:  # Edge is within the neighbor set
                        local_w = vertex_map[w]
                        subgraph_adj[local_u].append(local_w)
            
            # Find k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity._find_k_cores(subgraph_adj, k)
            sd[v] = len(k_cores)
        
        return sd
    
    @staticmethod
    def _find_k_cores(adj_list, k):
        """
        Find all k-cores in a graph represented by adjacency list.
        Returns a list of k-cores, where each k-core is a set of vertex indices.
        """
        n = len(adj_list)
        if n == 0:
            return []
        if k < 0:
            return []
        
        # Special case for k=0: find all connected components
        if k == 0:
            visited = [False] * n
            components = []
            for start in range(n):
                if visited[start]:
                    continue
                component = []
                queue = deque([start])
                visited[start] = True
                while queue:
                    v = queue.popleft()
                    component.append(v)
                    for u in adj_list[v]:
                        if not visited[u]:
                            visited[u] = True
                            queue.append(u)
                if component:
                    components.append(set(component))
            return components
        
        # For k > 0: compute k-core decomposition
        degree = [len(adj_list[v]) for v in range(n)]
        removed = [False] * n
        
        # Remove vertices with degree < k iteratively
        queue = deque()
        for v in range(n):
            if degree[v] < k:
                queue.append(v)
                removed[v] = True
        
        while queue:
            v = queue.popleft()
            for u in adj_list[v]:
                if not removed[u]:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        removed[u] = True
        
        # Remaining vertices form the k-core
        remaining = [v for v in range(n) if not removed[v]]
        if not remaining:
            return []
        
        # Find connected components in the k-core
        visited = [False] * n
        k_cores = []
        
        for start in remaining:
            if visited[start]:
                continue
            
            # BFS to find connected component
            component = []
            queue = deque([start])
            visited[start] = True
            
            while queue:
                v = queue.popleft()
                component.append(v)
                for u in adj_list[v]:
                    if not removed[u] and not visited[u]:
                        visited[u] = True
                        queue.append(u)
            
            if component:
                k_cores.append(set(component))
        
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
