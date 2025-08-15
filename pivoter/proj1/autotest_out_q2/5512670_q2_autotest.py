#!/usr/bin/env python3
# Auto-generated for 5512670

STUDENT_ID = "5512670"
STUDENT_NAME = "Weiyuan Wang"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        
        for v in range(n):
            neighbors = set(G.adj_list[v])
            if not neighbors:
                sd[v] = 0
                continue
                
            subgraph_nodes = neighbors
            subgraph_adj = {}
            node_degree = {}
            
            for u in subgraph_nodes:
                subgraph_adj[u] = []
                node_degree[u] = 0
            
            for u in subgraph_nodes:
                for neighbor in G.adj_list[u]:
                    if neighbor in subgraph_nodes:
                        subgraph_adj[u].append(neighbor)
                        node_degree[u] += 1
            
            core_numbers = kCoreBaseStructuralDiversity.compute_k_cores(subgraph_adj, node_degree)
            
            visited = set()
            k_core_count = 0
            
            for node in subgraph_nodes:
                if node not in visited and core_numbers.get(node, 0) >= k:
                    queue = deque()
                    queue.append(node)
                    visited.add(node)
                    is_k_core_component = True
                    
                    while queue:
                        current = queue.popleft()
                        if core_numbers.get(current, 0) < k:
                            is_k_core_component = False
                        
                        for neighbor in subgraph_adj[current]:
                            if neighbor not in visited and core_numbers.get(neighbor, 0) >= k:
                                visited.add(neighbor)
                                queue.append(neighbor)
                    
                    if is_k_core_component:
                        k_core_count += 1
            
            sd[v] = k_core_count
        
        return sd
    
    @staticmethod
    def compute_k_cores(adj_list, node_degree):
        """
        Compute core numbers for all nodes using the Batagelj-Zaversnik algorithm
        """
        degrees = node_degree.copy()
        nodes = list(degrees.keys())
        max_degree = max(degrees.values()) if degrees else 0
        
        core_numbers = {node: 0 for node in nodes}
        
        bins = [[] for _ in range(max_degree + 1)]
        for node, degree in degrees.items():
            bins[degree].append(node)
        
        for current_degree in range(max_degree + 1):
            while bins[current_degree]:
                node = bins[current_degree].pop()
                core_numbers[node] = current_degree
                
                for neighbor in adj_list[node]:
                    if degrees[neighbor] > current_degree:
                        index = degrees[neighbor]
                        bins[index].remove(neighbor)
                        degrees[neighbor] -= 1
                        bins[degrees[neighbor]].append(neighbor)
        
        return core_numbers

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
