#!/usr/bin/env python3
# Auto-generated for 5460237

STUDENT_ID = "5460237"
STUDENT_NAME = "Yuyang Sun"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        
        # Preprocess the neighbor sets of all vertices
        neighbors = [set(G.adj_list[v]) for v in range(n)]
        
        # Traverse all vertices and calculate the k-core based structural diversity of each vertex
        for v in range(n):
            deg_v = len(G.adj_list[v])
            if deg_v < k:
                sd[v] = 0
                continue
            
            # Constructing neighbor-induced subgraphs
            neighbor_list = list(neighbors[v])
            neighbor_set = neighbors[v]
            
            # Constructing adjacency lists of neighbor-induced subgraphs
            sub_adj = {
                u: [w for w in G.adj_list[u] if w in neighbor_set]
                for u in neighbor_list
            }

            
            # Calculate the number of k-cores in the neighbor-induced subgraph
            sd[v] = kCoreBaseStructuralDiversity._count_k_cores(sub_adj, neighbor_list, k)
        
        return sd

    @staticmethod
    def _count_k_cores(sub_adj, vertices, k):
        """
        Calculate the number of connected components of the k-core in the subgraph

        Parameters:
        sub_adj: adjacency list of the subgraph
        vertices: list of vertices in the subgraph
        k: parameter of the k-core
        
        Returns:
        int: number of connected components of the k-core
        """
        if not vertices or k <= 0:
            return 0
        
        # Calculate the degree of each vertex
        degrees = {u: len(sub_adj[u]) for u in vertices}
        
        # k-core decomposition: remove vertices with degree less than k
        remaining = set(vertices)
        queue = deque()
        
        # Initialization: Find all vertices with degree less than k
        for u in vertices:
            if degrees[u] < k:
                queue.append(u)
        
        # Iteratively remove vertices with degree less than k
        while queue:
            u = queue.popleft()
            if u not in remaining:
                continue
                
            remaining.remove(u)
            
            # Update neighbor's degree
            for v in sub_adj[u]:
                if v in remaining:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)
        
        # If no vertices remain, there is no k-core
        if not remaining:
            return 0
        
        # Calculate the number of connected components in the k-core
        visited = set()
        component_count = 0
        
        for u in remaining:
            if u not in visited:
                component_count += 1
                # BFS traverses connected components
                bfs_queue = deque([u])
                visited.add(u)
                
                while bfs_queue:
                    curr = bfs_queue.popleft()
                    for neighbor in sub_adj[curr]:
                        if neighbor in remaining and neighbor not in visited:
                            visited.add(neighbor)
                            bfs_queue.append(neighbor)
        
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
