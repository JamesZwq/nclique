#!/usr/bin/env python3
# Auto-generated for 5399710

STUDENT_ID = "5399710"
STUDENT_NAME = "Qinan Ji"

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

        # Iterate through all vertices and compute the number of connected k-core components
        # in each vertex's neighbor-induced subgraph
        for v in range(n):
            neighbours = G.adj_list[v]
            if not neighbours:
                sd[v] = 0
                continue

            # Build neighbor-induced subgraph (excluding vertex v itself)
            sub_G = {u: [] for u in neighbours}
            for u in neighbours:
                for neighbour_of_u in G.adj_list[u]:
                    # Only consider vertices in the neighbors set
                    if neighbour_of_u in sub_G:
                        sub_G[u].append(neighbour_of_u)

            # Compute the number of connected k-core components in the subgraph
            sd[v] = kCoreBaseStructuralDiversity._compute_k_core(sub_G, k)

        return sd

    @staticmethod
    def _compute_k_core(sub_G, k):
        """
        Compute the number of connected k-core components in the given subgraph
        """
        if not sub_G:
            return 0

        # If k <= 0, each vertex forms a 0-core, but only connected components count
        if k <= 0:
            return kCoreBaseStructuralDiversity._count_connected_components(sub_G)

        # First find the k-core (remove vertices with degree less than k)
        degrees = {u: len(neighbours) for u, neighbours in sub_G.items()}
        remaining = set(sub_G.keys())

        # Use queue to process vertices with degree less than k
        queue = deque()
        for u in sub_G:
            if degrees[u] < k:
                queue.append(u)

        # Iteratively remove vertices with degree less than k
        while queue:
            u = queue.popleft()
            if u not in remaining:
                continue

            remaining.remove(u)

            # Update degrees of neighbors
            for v in sub_G[u]:
                if v in remaining:
                    degrees[v] -= 1
                    if degrees[v] < k and v not in queue:
                        queue.append(v)

        # Build k-core subgraph
        k_core_subgraph = {}
        for u in remaining:
            k_core_subgraph[u] = [v for v in sub_G[u] if v in remaining]

        return kCoreBaseStructuralDiversity._count_connected_components(k_core_subgraph)

    @staticmethod
    def _count_connected_components(graph):
        """
        Compute the number of connected components in the graph using BFS.
        """
        # Handle empty graph case
        if not graph:
            return 0

        visited = set()
        components = 0

        # Iterate through all vertices
        for start_vertex in graph:

            # If not visited, it belongs to a new component
            if start_vertex not in visited:
                # Use BFS to traverse
                queue = deque([start_vertex])
                visited.add(start_vertex)

                # Continue BFS until all vertices in this component are visited
                while queue:
                    vertex = queue.popleft()

                    # Examine all neighbors of the current vertex
                    for neighbor in graph[vertex]:
                        # If neighbor hasn't been visited, add it to the current component
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

                components += 1

        return components


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
