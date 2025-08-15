#!/usr/bin/env python3
# Auto-generated for 5465844

STUDENT_ID = "5465844"
STUDENT_NAME = "Yinuo Yang"

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
        # if some vertex has degree less than k, then its structural diversity is 0
        for v in range(n):
            # beacause the graph have some duplicated edges, we need to use set to remove duplicates
            neighbors = set(G.adj_list[v])
            if not neighbors:
                sd[v] = 0
                continue
            sub_graph = {graph: [] for graph in neighbors}
            # find the subgraph induced by neighbors of v
            for neighbor in neighbors:
                for u in G.adj_list[neighbor]:
                    if u in neighbors and u != v:
                        sub_graph[neighbor].append(u)
            if not sub_graph:
                return 0
            sd[v] = kCoreBaseStructuralDiversity.compute_k_core_diversity(sub_graph, k)

        return sd
        
    @staticmethod
    def compute_k_core_diversity(sub_graph, k):        
        degree = {graph: len(neighbors) for graph, neighbors in sub_graph.items()}
        remove = set()

        queue = deque()
        # Initialize the queue with vertices of degree less than k
        for graph in sub_graph:
            if degree[graph] < k:
                queue.append(graph)
                remove.add(graph)
        
        # Process the queue, removing vertices and updating degrees
        while queue:
            current = queue.popleft()
            for neighbor in sub_graph[current]:
                if neighbor not in remove:
                    degree[neighbor] -= 1
                    if degree[neighbor] < k:
                        queue.append(neighbor)
                        remove.add(neighbor)
        # The remaining vertices in the subgraph are those with degree >= k
        remaining_vertices = [v for v in sub_graph if v not in remove]

        if not remaining_vertices:
            return 0
        # Count the number of connected components in the remaining subgraph
        num_k_core = kCoreBaseStructuralDiversity.dfs(remaining_vertices,remove,sub_graph)
        return num_k_core
    
    @staticmethod
    def dfs(remaining_vertices, remove, sub_graph):
        # Perform DFS to count connected components in the subgraph that are not removed
        visited = set()
        num_k_core = 0
        stack = []        
        for vertex in remaining_vertices:
            if vertex not in visited:
                num_k_core += 1
                visited.add(vertex)
                stack.append(vertex)
                while stack:
                    current = stack.pop()
                    for neighbor in sub_graph[current]:
                        if neighbor not in visited and neighbor not in remove:
                            visited.add(neighbor)
                            stack.append(neighbor)    
        return num_k_core


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
