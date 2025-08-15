#!/usr/bin/env python3
# Auto-generated for 5486847

STUDENT_ID = "5486847"
STUDENT_NAME = "Mengda Wang"

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
            # Get the neighborhood of vertex v
            neighbors = G.adj_list[v]
            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # Build the neighborhood-induced subgraph
            neighbor_set = set(neighbors)

            # Calculate the number of k-cores in the neighborhood-induced subgraph
            sd[v] = kCoreBaseStructuralDiversity._count_connected_k_cores(G, neighbor_set, k)

        return sd

    @staticmethod
    def _count_connected_k_cores(G, neighbor_set, k):
        """
        Calculate the number of connected k-cores in the neighborhood-induced subgraph
        """
        if len(neighbor_set) == 0:
            return 0

        # Build adjacency list for the neighborhood-induced subgraph
        subgraph_adj = {}
        for v in neighbor_set:
            subgraph_adj[v] = []

        # Extract edges from the original graph for the neighborhood subgraph
        for v in neighbor_set:
            for neighbor in G.adj_list[v]:
                if neighbor in neighbor_set:
                    subgraph_adj[v].append(neighbor)

        # Calculate the number of connected k-cores
        return kCoreBaseStructuralDiversity._find_connected_k_cores_count(subgraph_adj, k)

    @staticmethod
    def _find_connected_k_cores_count(adj_list, k):
        """
        Find all connected k-cores in the graph and return the count
        """
        if not adj_list:
            return 0

        # Copy the adjacency list
        graph = {v: list(neighbors) for v, neighbors in adj_list.items()}

        # First perform k-core decomposition
        k_core_vertices = kCoreBaseStructuralDiversity._extract_k_core_vertices(graph, k)

        if not k_core_vertices:
            return 0

        # Find connected components in the k-core vertices
        return kCoreBaseStructuralDiversity._count_connected_components(graph, k_core_vertices)

    @staticmethod
    def _extract_k_core_vertices(graph, k):
        """
        Extract the vertex set of k-core from the graph
        """
        if not graph:
            return set()

        # Use queue for k-core decomposition
        queue = deque()
        degrees = {}

        # Calculate initial degrees
        for v in graph:
            degrees[v] = len(graph[v])
            if degrees[v] < k:
                queue.append(v)

        # Remove vertices with degree less than k
        while queue:
            v = queue.popleft()
            if v not in graph:
                continue

            # Remove vertex v
            for neighbor in graph[v]:
                if neighbor in graph:
                    graph[neighbor].remove(v)
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k and neighbor not in queue:
                        queue.append(neighbor)

            del graph[v]
            del degrees[v]

        # Return the remaining k-core vertices
        return set(graph.keys())

    @staticmethod
    def _count_connected_components(graph, vertices):
        """
        Count the number of connected components in the given vertex set
        """
        if not vertices:
            return 0

        visited = set()
        components = 0

        for v in vertices:
            if v not in visited:
                # Use BFS to find connected component
                queue = deque([v])
                visited.add(v)

                while queue:
                    current = queue.popleft()
                    for neighbor in graph[current]:
                        if neighbor in vertices and neighbor not in visited:
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
