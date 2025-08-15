#!/usr/bin/env python3
# Auto-generated for 5452768

STUDENT_ID = "5452768"
STUDENT_NAME = "Zhiyuan Yan"

# ======= 学生代码 =======
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The input graph
        k : int
            The core number parameter
            
        Returns
        -------
        List[int]
            τ_k(v) for all vertices v in G
        """
        n = G.vertex_num
        τ = [0] * n
        
        # Process each vertex
        for v in range(n):
            if not G.adj_list[v]:
                continue
                
            # Get neighbor set N(v)
            neighbors = set(G.adj_list[v])
            
            # Build neighbor-induced subgraph G[N(v)]
            # Create adjacency list for the induced subgraph
            neighbor_map = {u: i for i, u in enumerate(neighbors)}
            sub_n = len(neighbors)
            sub_adj = [[] for _ in range(sub_n)]
            
            # Add edges within the neighbor set
            for u in neighbors:
                u_idx = neighbor_map[u]
                for w in G.adj_list[u]:
                    if w in neighbors and w != u:
                        w_idx = neighbor_map[w]
                        sub_adj[u_idx].append(w_idx)
            
            # Remove duplicates and self-loops
            for i in range(sub_n):
                sub_adj[i] = list(set(sub_adj[i]))
            
            # Find k-cores in the neighbor-induced subgraph
            τ[v] = kCoreBaseStructuralDiversity._count_k_cores(sub_adj, k)
        
        return τ
    
    @staticmethod
    def _count_k_cores(adj_list, k):
        n = len(adj_list)
        if n == 0:
            return 0
        
        degree = [len(adj_list[i]) for i in range(n)]
        
        # Find vertices in k-core using iterative pruning
        in_kcore = [True] * n
        queue = []
        
        # Initial pruning: remove vertices with degree < k
        for v in range(n):
            if degree[v] < k:
                queue.append(v)
                in_kcore[v] = False
        
        # Iteratively remove vertices with degree < k
        while queue:
            v = queue.pop(0)
            for u in adj_list[v]:
                if in_kcore[u]:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        in_kcore[u] = False
        
        # Count connected components in the k-core
        visited = [False] * n
        num_components = 0
        
        def dfs(v):
            visited[v] = True
            for u in adj_list[v]:
                if in_kcore[u] and not visited[u]:
                    dfs(u)
        
        # Count connected components
        for v in range(n):
            if in_kcore[v] and not visited[v]:
                dfs(v)
                num_components += 1
        
        return num_components

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
