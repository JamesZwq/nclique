#!/usr/bin/env python3
# Auto-generated for 5524919

STUDENT_ID = "5524919"
STUDENT_NAME = "Yulin Song"

# ======= 学生代码 =======
from collections import deque

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

        Traverse all vertices and compute the number of k-cores in the neighbor-induced subgraph of each vertex.
        """
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            # Obtain the neighbors of the vertex and construct the corresponding subgraph.
            neighbours = set(G.adj_list[v])
            if not neighbours:
                sd[v] = 0
                continue
            # Note that the subgraph does not include the vertex v itself; we construct it using an adjacency list.
            sub_G = kCoreBaseStructuralDiversity._build_neighbor_induced_subgraph(G, neighbours)
            # use compute k core to count sub graph's k-core num
            sd[v] = kCoreBaseStructuralDiversity._count_k_core(sub_G, k)


        return sd

    @staticmethod
    def _build_neighbor_induced_subgraph(G, neighbours):
        """
        Construct the adjacency list of the neighbor-induced subgraph (excluding the vertex itself).
        Parameters:
	      •	G: the original graph
	      •	neighbours: the set of neighbors of vertex v
        Returns:
	    •	dict: the adjacency list of the induced subgraph
        """
        sub_G = {u: [] for u in neighbours}
        for u in neighbours:
            for neighbours_of_u in G.adj_list[u]:
                if neighbours_of_u in neighbours:
                    sub_G[u].append(neighbours_of_u)
        return sub_G

    @staticmethod
    def _count_k_core(sub_G, k):
        """
        Compute the number of connected components within the k-core of the subgraph.
        Since the value of k is known, we simply compare each vertex’s degree with k,
        remove vertices with degrees less than k, and update the degrees of their neighbors accordingly.
        Initially, all vertices with degree < k are collected and processed together.
        During this process, some neighbors’ degrees may also drop below k, and these vertices must be processed as well.
        The above process is the k-core computation step.
        Next is the counting step:
        If the remaining vertices (remains) are not empty, the k-core exists. We then count how many connected components are
        present in the remaining subgraph using BFS or DFS.
        Only vertices in the k-core are included in the connectivity check.

        sub_G: dict, adjacency list of the subgraph
        k: int, the k value for the k-core

        int, the number of connected components in the k-core
        """
        if not sub_G:
            return 0
        degrees = {u: len(neighbours) for u, neighbours in sub_G.items()}
        deleted = set()
        vertices = deque()
        # First, identify and process all vertices with a degree less than k.
        for u in sub_G:
            if degrees[u] < k:
                vertices.append(u)
                deleted.add(u)
        # Process these vertices, and during this process, some other vertices may have their degrees reduced below k.
        while vertices:
            u = vertices.popleft()
            for v in sub_G[u]:
                if v not in deleted:
                    degrees[v] -= 1
                    # Check whether the degree is less than k; if so, add the vertex to the processing queue.
                    if degrees[v] < k:
                        deleted.add(v)
                        vertices.append(v)
        remains = [u for u in sub_G if u not in deleted]
        if not remains:
            return 0
        # If the remains set is not empty, then a k-core exists. Count how many k-cores there are using BFS.
        visited = set()
        k_core_count = 0
        for node in remains:
            if node not in visited:
                k_core_count += 1
                queue = deque()
                queue.append(node)
                visited.add(node)
                while queue:
                    u = queue.popleft()
                    for v in sub_G[u]:
                        # Perform connectivity checks only on the vertices within the k-core.
                        if v not in deleted and v not in visited:
                            visited.add(v)
                            queue.append(v)
        return k_core_count

    ################################################################################
    # You can only define any auxiliary functions in class

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
