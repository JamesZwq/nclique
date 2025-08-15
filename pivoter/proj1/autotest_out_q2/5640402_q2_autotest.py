#!/usr/bin/env python3
# Auto-generated for 5640402

STUDENT_ID = "5640402"
STUDENT_NAME = "Hao Zhang"

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

        # Iterate over all vertices
        for v in range(n):
            # Extract neighbors of node v
            neighbors = G.adj_list[v]
            if not neighbors:
                sd[v] = 0
                continue

            # Build the neighbor-induced subgraph (excluding v itself)
            neighbors = set(neighbors)  # Ensure uniqueness

            # Build the neighbor-induced subgraph (adjacency list)
            subgraph_G = {u: [] for u in neighbors}
            for u in neighbors:
                for neighbor_of_u in G.adj_list[u]:
                    if neighbor_of_u in neighbors:
                        subgraph_G[u].append(neighbor_of_u)

            # Compute the number of connected components in the k-core of the subgraph
            sd[v] = kCoreBaseStructuralDiversity._compute_k_core(subgraph_G, k)

        return sd

    @staticmethod
    def _compute_k_core(subgraph_G, k):
        """
        Compute the k-core of a given subgraph.
        """
        if not subgraph_G:
            return 0

        degree = {u: len(neighbors) for u, neighbors in subgraph_G.items()}
        deleted = set()

        # Prune all nodes with degree < k
        vertices = deque()
        for u in subgraph_G:
            if degree[u] < k:
                vertices.append(u)
                deleted.add(u)

        # Iteratively prune neighbors whose degree falls below k
        while vertices:
            u = vertices.popleft()
            for v in subgraph_G[u]:
                if v not in deleted:
                    degree[v] -= 1
                    # If the degree of v is less than k, add it to the list of deleted
                    if degree[v] < k:
                        vertices.append(v)
                        deleted.add(v)

        # Collect remaining nodes after pruning
        remaining = [u for u in subgraph_G if u not in deleted]
        if not remaining:
            return 0

        # Use BFS to count connected components in the remaining graphe
        visited = set()
        k_core_count = 0

        for node in remaining:
            if node not in visited:
                k_core_count += 1
                queue = deque()
                queue.append(node)
                visited.add(node)

                while queue:
                    u = queue.popleft()
                    for v in subgraph_G[u]:
                        # Only visit non-deleted, unvisited nodes
                        if v not in deleted and v not in visited:
                            visited.add(v)
                            queue.append(v)

        return k_core_count



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
