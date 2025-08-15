#!/usr/bin/env python3
# Auto-generated for 5129375

STUDENT_ID = "5129375"
STUDENT_NAME = "Peiying Tan"

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
            # Get neighbors of v
            neighbors = G.adj_list[v]
            if not neighbors:
                sd[v] = 0
                continue

            # Create neighbor-induced subgraph
            subgraph_nodes = set(neighbors)
            subgraph_adj = {}
            degree = {}

            # Build adjacency list and degree for the subgraph
            for u in neighbors:
                subgraph_adj[u] = []
                degree[u] = 0

            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in subgraph_nodes:
                        subgraph_adj[u].append(w)
                        degree[u] += 1

            # Compute k-core decomposition of the neighbor-induced subgraph
            # Using the standard k-core algorithm (peeling process)
            queue = deque()
            remaining_nodes = set(subgraph_nodes)

            # Initialize queue with nodes having degree < k
            for u in subgraph_nodes:
                if degree[u] < k:
                    queue.append(u)

            # Peeling process
            while queue:
                u = queue.popleft()
                if u not in remaining_nodes:
                    continue

                remaining_nodes.remove(u)

                for w in subgraph_adj[u]:
                    if w in remaining_nodes:
                        degree[w] -= 1
                        if degree[w] < k:
                            queue.append(w)

            # Now count the number of connected components in the remaining graph (k-cores)
            visited = set()
            component_count = 0

            for u in remaining_nodes:
                if u not in visited:
                    # BFS to find connected component
                    component_count += 1
                    q = deque([u])
                    visited.add(u)

                    while q:
                        current = q.popleft()
                        for w in subgraph_adj[current]:
                            if w in remaining_nodes and w not in visited:
                                visited.add(w)
                                q.append(w)

            sd[v] = component_count
        return sd


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
