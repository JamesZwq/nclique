#!/usr/bin/env python3
# Auto-generated for 5513586

STUDENT_ID = "5513586"
STUDENT_NAME = "Shixu Yang"

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
        # Traverse all vertices and calculate the number of k-cores in the neighbor-induced-subgraph of each vertex
        n = G.vertex_num
        if n == 0:
            return []

        # Directly return a list of τ values for each vertex
        tau_values = []
        for v in range(n):
            tau_v = kCoreBaseStructuralDiversity._compute_tau_for_vertex(G, v, k)
            tau_values.append(tau_v)

        return tau_values

    @staticmethod
    def _compute_tau_for_vertex(G, vertex, k):
        neighbors = set(G.adj_list[vertex])

        if len(neighbors) < k:
            return 0  # Insufficient number of neighbors to form a k-core

        # Neighbor induced subgraph G[N(v)]
        neighbor_subgraph = kCoreBaseStructuralDiversity._build_neighbor_induced_subgraph(G, neighbors)

        # Call compute k-core to calculate the number of k-cores in the sub graph
        k_core_vertices = kCoreBaseStructuralDiversity._compute_k_core_vertices(neighbor_subgraph, k)

        if not k_core_vertices:
            return 0

        # If the remaining vertices are not empty, there is a k-core. Use BFS to query how many k-cores there are.
        component_count = kCoreBaseStructuralDiversity._count_connected_components_in_k_core(neighbor_subgraph, k_core_vertices)

        return component_count

    @staticmethod
    def _build_neighbor_induced_subgraph(G, neighbors):
        # Construct neighbor induced subgraph G[N(v)]
        # Create a mapping of neighbors to new IDs
        neighbor_list = list(neighbors)
        neighbor_to_id = {neighbor: i for i, neighbor in enumerate(neighbor_list)}

        # Constructing the adjacency list of the subgraph
        subgraph_adj = [[] for _ in range(len(neighbors))]

        # Adding edges between neighbors
        for neighbor in neighbors:
            for adj_neighbor in G.adj_list[neighbor]:
                if adj_neighbor in neighbors and adj_neighbor != neighbor:
                    # Avoid duplicate edges
                    if neighbor_to_id[adj_neighbor] not in subgraph_adj[neighbor_to_id[neighbor]]:
                        subgraph_adj[neighbor_to_id[neighbor]].append(neighbor_to_id[adj_neighbor])
                        subgraph_adj[neighbor_to_id[adj_neighbor]].append(neighbor_to_id[neighbor])

        return subgraph_adj

    @staticmethod
    def _compute_k_core_vertices(subgraph_adj, k):
        """
        Find all vertices with degree < k
        Update the degree of other vertices
        """
        n = len(subgraph_adj)
        if n == 0:
            return []

        # Calculate the degree of each vertex
        degrees = [len(adj_list) for adj_list in subgraph_adj]
        in_k_core = [True] * n

        # Use a queue to store vertices with degree less than k
        queue = deque()
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)
                in_k_core[i] = False

        # Iteratively remove vertices with degree less than k
        while queue:
            u = queue.popleft()
            for v in subgraph_adj[u]:
                if in_k_core[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)
                        in_k_core[v] = False

        # Returns the vertices in the k-core
        return [i for i in range(n) if in_k_core[i]]

    @staticmethod
    def _count_connected_components_in_k_core(subgraph_adj, k_core_vertices):
        # Use BFS to calculate the number of connected components in k-core
        if not k_core_vertices:
            return 0

        n = len(subgraph_adj)
        visited = [False] * n
        component_count = 0

        # Convert k-core vertices into sets
        k_core_set = set(k_core_vertices)

        for start_vertex in k_core_vertices:
            if not visited[start_vertex]:
                # Use BFS to traverse connected components
                queue = deque([start_vertex])
                visited[start_vertex] = True

                while queue:
                    u = queue.popleft()
                    for v in subgraph_adj[u]:
                        if v in k_core_set and not visited[v]:
                            visited[v] = True
                            queue.append(v)

                component_count += 1

        return component_count

    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################


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
