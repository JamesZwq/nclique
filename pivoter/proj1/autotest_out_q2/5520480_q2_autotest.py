#!/usr/bin/env python3
# Auto-generated for 5520480

STUDENT_ID = "5520480"
STUDENT_NAME = "Sridhar Surya - Revathi Sridhar Surya"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules
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
        tau = [0] * n

        for v in range(n):
            # Get neighbors of vertex v
            neighbors = G.adj_list[v]

            if len(neighbors) == 0:
                # Isolated vertex has no k-cores in its neighbor-induced subgraph
                tau[v] = 0
                continue

            # Create neighbor-induced subgraph
            # Map original vertex IDs to local IDs (0, 1, 2, ...)
            neighbor_to_local = {neighbor: i for i, neighbor in enumerate(neighbors)}
            local_adj_list = [[] for _ in range(len(neighbors))]

            # Build adjacency list for neighbor-induced subgraph
            for i, u in enumerate(neighbors):
                for neighbor_of_u in G.adj_list[u]:
                    if neighbor_of_u in neighbor_to_local:
                        local_v = neighbor_to_local[neighbor_of_u]
                        local_adj_list[i].append(local_v)

            # Find k-cores in the neighbor-induced subgraph
            tau[v] = kCoreBaseStructuralDiversity._count_k_cores(local_adj_list, k)

        return tau

    @staticmethod
    def _count_k_cores(adj_list, k):
        """
        Count the number of k-cores in a graph represented by adj_list
        """
        n = len(adj_list)
        if n == 0:
            return 0

        # Step 1: Remove vertices with degree < k iteratively until convergence
        active = [True] * n
        degrees = [len(adj_list[i]) for i in range(n)]
        queue = deque()

        # Initialize queue with vertices having degree < k
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)
                active[i] = False

        # Remove vertices with insufficient degree
        while queue:
            u = queue.popleft()
            for v in adj_list[u]:
                if active[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)
                        active[v] = False

        # Step 2: Find connected components in the remaining graph
        visited = [False] * n
        k_core_count = 0

        for i in range(n):
            if active[i] and not visited[i]:
                # Found a new k-core component
                k_core_count += 1
                # BFS to mark all vertices in this component
                bfs_queue = deque([i])
                visited[i] = True

                while bfs_queue:
                    u = bfs_queue.popleft()
                    for v in adj_list[u]:
                        if active[v] and not visited[v]:
                            visited[v] = True
                            bfs_queue.append(v)

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
