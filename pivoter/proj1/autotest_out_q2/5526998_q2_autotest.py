#!/usr/bin/env python3
# Auto-generated for 5526998

STUDENT_ID = "5526998"
STUDENT_NAME = "Jingquan Mao"

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
        """
        # TODO
        n = G.vertex_num
        structural_diversity = [0] * n

        for v in range(n):
            # Get neighbors of vertex v
            neighbors = G.adj_list[v]
            if not neighbors:
                continue

            # Build neighbor-induced subgraph
            neighbor_set = set(neighbors)
            neighbor_subgraph = {}
            for u in neighbors:
                # Get neighbors of u that are also in v's neighbor set
                u_neighbors = [w for w in G.adj_list[u] if w in neighbor_set]
                neighbor_subgraph[u] = u_neighbors

            # Calculate k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity._compute_k_cores(neighbor_subgraph, k)
            structural_diversity[v] = len(k_cores)

        return structural_diversity

    @staticmethod
    def _compute_k_cores(graph, k):
        """
        Compute all k-cores in a graph and return their connected components.

        Args:
            graph: Dictionary representing the adjacency list of the input graph
            k: Integer specifying the k value for k-core calculation

        Returns:
            A list of lists, where each inner list contains vertices of a k-core connected component
        """
        if not graph:
            return []

        # Create a copy of the graph and calculate initial degrees
        subgraph = {u: set(neighbors) for u, neighbors in graph.items()}
        degrees = {u: len(neighbors) for u, neighbors in subgraph.items()}

        # Peeling algorithm: remove vertices with degree < k iteratively
        queue = deque([u for u in degrees if degrees[u] < k])
        removed = set()

        while queue:
            u = queue.popleft()
            if u in removed:
                continue
            removed.add(u)

            # Update degrees of neighbors
            for v in subgraph[u]:
                if v not in removed:
                    degrees[v] -= 1
                    if degrees[v] < k and v not in queue:
                        queue.append(v)

        # Remaining vertices form the candidate set for k-cores
        remaining = [u for u in graph if u not in removed]
        if not remaining:
            return []

        # Find all connected components in the remaining subgraph
        visited = set()
        components = []

        for u in remaining:
            if u not in visited:
                # BFS to find connected component
                component = []
                q = deque([u])
                visited.add(u)

                while q:
                    current = q.popleft()
                    component.append(current)

                    for v in subgraph[current]:
                        if v in remaining and v not in visited:
                            visited.add(v)
                            q.append(v)

                components.append(component)

        return components

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
