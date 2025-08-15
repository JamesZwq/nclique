#!/usr/bin/env python3
# Auto-generated for 5504341

STUDENT_ID = "5504341"
STUDENT_NAME = "Jiahang Zhang"

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
            
            neighbors = list(set(G.adj_list[v]))  
            if not neighbors:  
                sd[v] = 0
                continue
            
            # 构建邻居诱导子图的邻接表
            neighbor_adj = kCoreBaseStructuralDiversity._build_neighbor_adj(G, neighbors)
            
            # 计算邻居诱导子图的k-core
            k_core = kCoreBaseStructuralDiversity._compute_k_core(neighbor_adj, k)
            
            # 计算k-core中的连通分量数量
            sd[v] = kCoreBaseStructuralDiversity._count_connected_components(k_core, neighbor_adj)
            
        return sd

    @staticmethod
    def _build_neighbor_adj(G, neighbors):
        neighbor_set = set(neighbors)
        neighbor_adj = {u: [] for u in neighbors}
        
        for u in neighbors:
            
            for v in G.adj_list[u]:
                if v in neighbor_set and v != u and v not in neighbor_adj[u]:
                    neighbor_adj[u].append(v)
                    neighbor_adj[v].append(u)  
        
        return neighbor_adj

    @staticmethod
    def _compute_k_core(adj, k):
        
        if not adj:
            return set()
            
        
        current_adj = {u: list(vs) for u, vs in adj.items()}
        degrees = {u: len(vs) for u, vs in current_adj.items()}
        in_core = {u: True for u in current_adj}
        queue = deque([u for u in current_adj if degrees[u] < k])
        
        while queue:
            u = queue.popleft()
            if not in_core[u]:  
                continue
            in_core[u] = False  
            
            for v in current_adj[u]:
                if in_core[v]:
                    degrees[v] -= 1
                    if degrees[v] < k and v not in queue:
                        queue.append(v)
        
        return {u for u in in_core if in_core[u]}

    @staticmethod
    def _count_connected_components(k_core, adj):
        if not k_core:
            return 0
            
        visited = set()
        count = 0
        core_adj = {u: [v for v in adj[u] if v in k_core] for u in k_core}
        
        for u in k_core:
            if u not in visited:
                count += 1
                queue = deque([u])
                visited.add(u)
                
                while queue:
                    current = queue.popleft()
                    for neighbor in core_adj[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
        
        return count

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
