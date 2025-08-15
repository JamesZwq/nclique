#!/usr/bin/env python3
# Auto-generated for 5602135

STUDENT_ID = "5602135"
STUDENT_NAME = "Yihong Zhang"

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

        # Process vertices
        for v in range(n):
            # Get neighbors of vertex v
            neighbors = set(G.adj_list[v])

            if not neighbors:
                # No neighbors, so no k-cores
                sd[v] = 0
                continue

            # Create neighbor-induced subgraph
            neighbor_list = list(neighbors)
            neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbor_list)}

            # Build adjacency list for neighbor-induced subgraph
            neighbor_adj = [[] for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list):
                for neighbor in G.adj_list[u]:
                    # Only include edges within neighbors
                    if neighbor in neighbors:
                        neighbor_adj[i].append(neighbor_to_idx[neighbor])

            # Find all k-cores in the neighbor-induced subgraph
            sd[v] = kCoreBaseStructuralDiversity._count_k_cores(neighbor_adj, k)

        return sd

    @staticmethod
    def _count_k_cores(adj_list, k):
        """
        Count the number of k-cores in the given graph represented by adj_list.
        :param adj_list: Adjacency list representation of the graph
        :param k: The k value for k-core
        :return: Number of k-cores in the graph
        """
        # Empty graph
        if not adj_list:
            return 0

        # Compute k-core decomposition
        degrees = [len(neighbors) for neighbors in adj_list]
        removed = [False] * len(adj_list)

        # Use a queue to process vertices with degree < k
        queue = deque()
        for i in range(len(adj_list)):
            if degrees[i] < k:
                # Push to queue
                queue.append(i)

        # Remove vertices with degree < k iteratively
        while queue:
            u = queue.popleft()
            # Already removed
            if removed[u]:
                continue
            removed[u] = True

            # Update degrees of neighbors
            for v in adj_list[u]:
                if removed[v]:
                    continue

                degrees[v] -= 1
                if degrees[v] < k:
                    queue.append(v)

        # Find connected components in the remaining graph (k-core)
        visited = [False] * len(adj_list)
        num_k_cores = 0

        for i in range(len(adj_list)):
            # Skip removed or visited
            if removed[i] or visited[i]:
                continue

            # Found a new k-core component
            num_k_cores += 1
            # BFS to mark all vertices in this component
            bfs_queue = deque([i])
            visited[i] = True

            while bfs_queue:
                u = bfs_queue.popleft()
                for v in adj_list[u]:
                    # Skip removed or visited
                    if removed[v] or visited[v]:
                        continue
                    visited[v] = True
                    # Add to BFS queue
                    bfs_queue.append(v)

        return num_k_cores

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
