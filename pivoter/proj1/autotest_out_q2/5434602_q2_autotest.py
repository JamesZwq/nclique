#!/usr/bin/env python3
# Auto-generated for 5434602

STUDENT_ID = "5434602"
STUDENT_NAME = "Mike Ling"

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
            The input graph represented as an adjacency list.
        k : int
            The minimum degree threshold for the k-core.

        Returns
        -------
        List[int]
            A list of integers where each value represents the k-core based
            structural diversity (τ_k(v)) for each vertex in the graph.
        """
        # Number of vertices in the graph
        n = G.vertex_num
        # Initialize the structural diversity array with zeros
        sd = [0] * n

        # Iterate over each vertex in the graph
        for v in range(n):
            # Get the neighbors of the current vertex
            neighbors = G.adj_list[v]
            # If the vertex has no neighbors, its structural diversity is 0
            if not neighbors:
                sd[v] = 0
                continue

            # Construct the ego network subgraph for the current vertex
            sub_G = {u: [] for u in neighbors}
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbors:
                        sub_G[u].append(w)

            # Compute the k-core size of the ego network
            sd[v] = kCoreBaseStructuralDiversity._k_core_size(sub_G, k)

        return sd

    @staticmethod
    def _k_core_size(sub_G, k):
        """
        Compute the size of the k-core for a given subgraph.

        Parameters
        ----------
        sub_G : dict
            The ego network subgraph represented as an adjacency list.
        k : int
            The minimum degree threshold for the k-core.

        Returns
        -------
        int
            The number of connected components in the k-core of the subgraph.
        """
        # If the subgraph is empty, return 0
        if not sub_G:
            return 0

        # Initialize the degree of each vertex in the subgraph
        degrees = {u: len(neighbors) for u, neighbors in sub_G.items()}
        # Set to keep track of deleted vertices
        deleted = set()

        # Queue to process vertices with degree less than k
        vertices = deque()
        for u in sub_G:
            if degrees[u] < k:
                vertices.append(u)
                deleted.add(u)

        # Iteratively remove vertices with degree less than k
        while vertices:
            u = vertices.popleft()
            for v in sub_G[u]:
                if v not in deleted:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        vertices.append(v)
                        deleted.add(v)

        # Get the remaining vertices in the k-core
        remains = [u for u in sub_G if u not in deleted]
        if not remains:
            return 0

        # Use BFS/DFS to count the number of connected components
        visited = set()
        k_core_size = 0

        for u in remains:
            if u not in visited:
                k_core_size += 1
                queue = deque()
                queue.append(u)
                visited.add(u)
                while queue:
                    v = queue.popleft()
                    for w in sub_G[v]:
                        if w not in deleted and w not in visited:
                            visited.add(w)
                            queue.append(w)

        return k_core_size

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
