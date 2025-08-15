#!/usr/bin/env python3
# Auto-generated for 5545575

STUDENT_ID = "5545575"
STUDENT_NAME = "Yuchen Shi"

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
        # Iterate through all vertices to calculate the k-core count in each vertex's neighbor-induced subgraph
        
        for v in range(n):
            # Get the neighbors of vertex v and form a subgraph

            neighbours = G.adj_list[v]
            if not neighbours:
                sd[v]= 0
                continue

            # Note: the subgraph does not include vertex v itself, we form an adjacency list to handle this
            neighbours = set(neighbours)  # Use set to remove duplicates
            sub_G = {u: [] for u in neighbours}
            for u in neighbours:
                for neighbours_of_u in G.adj_list[u]:
                    if neighbours_of_u in neighbours:
                        sub_G[u].append(neighbours_of_u)
            # Call compute k core to calculate the number of k-cores in the subgraph

            sd[v]= kCoreBaseStructuralDiversity._compute_k_core(sub_G, k)
            
        print()
        return sd

    @staticmethod
    def _compute_k_core(sub_G, k):
        """
        Calculate the number of connected components in the k-core of a subgraph
        
        Parameters
        ----------
        sub_G : dict
            Adjacency list representation of subgraph {vertex: [neighbors]}
        k : int
            k value
            
        Returns
        -------
        int : Number of connected components in the k-core
        """
        if not sub_G:
            return 0
            
        # If k=0, directly calculate the number of connected components
        if k == 0:
            return kCoreBaseStructuralDiversity._count_connected_components(sub_G)
        
        degrees = {v: len(neighbors) for v, neighbors in sub_G.items()}
        deleted = set()

        # Initialize queue with all vertices having degree less than k
        vertices = deque([v for v in sub_G if degrees[v] < k])
        
        while vertices:
            v = vertices.popleft()
            if v in deleted:
                continue
                
            deleted.add(v)
            # Update all neighbors of v
            for neighbor in sub_G[v]:
                if neighbor not in deleted:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        vertices.append(neighbor)

        remains= [v for v in sub_G if v not in deleted]

        if not remains:
            return 0
        
        visited = set()
        k_core_count = 0

        for node in remains:
            if node not in visited:
                k_core_count += 1
                # Depth-first search to traverse connected components
                queue=deque()
                queue.append(node)
                visited.add(node)
                while queue:
                    u = queue.popleft()
                    for v in sub_G[u]:
                        if v not in visited and v not in deleted:
                            visited.add(v)
                            queue.append(v)
        
        
        return k_core_count

    @staticmethod
    def _count_connected_components(graph):
        """
        Calculate the number of connected components in a graph
        
        Parameters
        ----------
        graph : dict
            Adjacency list representation of graph {vertex: [neighbors]}
            
        Returns
        -------
        int : Number of connected components
        """
        if not graph:
            return 0
            
        visited = set()
        component_count = 0
        
        for node in graph:
            if node not in visited:
                component_count += 1
                # Use BFS to traverse connected components
                queue = deque([node])
                visited.add(node)
                
                while queue:
                    current = queue.popleft()
                    for neighbor in graph[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
        
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
