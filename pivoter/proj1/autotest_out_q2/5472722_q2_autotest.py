#!/usr/bin/env python3
# Auto-generated for 5472722

STUDENT_ID = "5472722"
STUDENT_NAME = "Tao Zhou"

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
        List[int] # τ_k(v) for all v
        """
        n = G.vertex_num
        # Handle edge case: empty graph
        if n == 0:
            return []

        sd = []
        # Initialize diversity scores for all vertices
        for i in range(n):
            sd.append(0)

        # Calculate k-core structural diversity for each vertex
        v = 0
        while v < n:
            # Handle edge case: check if vertex exists in adjacency list
            if v >= len(G.adj_list):
                sd[v] = 0
                v += 1
                continue

            # Get neighbors of vertex v
            neighbors = G.adj_list[v]

            # Handle isolated vertices (no neighbors)
            if not neighbors:
                sd[v] = 0
                v += 1
                continue

            # Build neighbor induced subgraph
            neighbor_subgraph = kCoreBaseStructuralDiversity.build_neighbor_subgraph(G, neighbors)

            # Find k-core connected components in neighbor subgraph
            k_core_components = kCoreBaseStructuralDiversity.find_k_core_components(neighbor_subgraph, k)

            # Count the number of k-core connected components
            sd[v] = len(k_core_components)
            v += 1

        return sd

    @staticmethod
    def build_neighbor_subgraph(G, neighbors):
        """Build neighbor induced subgraph"""
        # Handle edge case: no neighbors
        if not neighbors:
            return {}

        # Put neighbors into dictionary for efficient lookup
        neighbor_dict = {}
        for neighbor in neighbors:
            neighbor_dict[neighbor] = True

        # Build adjacency list
        subgraph = {}
        for neighbor in neighbors:
            subgraph[neighbor] = []

        # Build adjacency relationships
        for neighbor in neighbors:
            # Handle edge case: check if neighbor exists in graph
            if neighbor < len(G.adj_list):
                for nn in G.adj_list[neighbor]:
                    if nn in neighbor_dict:
                        subgraph[neighbor].append(nn)

        return subgraph

    @staticmethod
    def find_k_core_components(subgraph, k):
        """Find k-core connected components in subgraph"""
        # Handle edge case: empty subgraph
        if not subgraph:
            return []

        # Handle edge case: k <= 0
        if k <= 0:
            # All nodes form one component if k <= 0
            if subgraph:
                return [set(subgraph.keys())]
            return []

        # First find k-core nodes
        k_core_nodes = kCoreBaseStructuralDiversity.k_core_decomposition(subgraph, k)

        if not k_core_nodes:
            return []

        # Build subgraph within k-core nodes
        k_core_subgraph = {}
        for node in k_core_nodes:
            k_core_subgraph[node] = []

        for node in k_core_nodes:
            for neighbor in subgraph[node]:
                if neighbor in k_core_nodes:
                    k_core_subgraph[node].append(neighbor)

        # Find connected components using BFS
        visited = set()
        components = []

        for node in k_core_nodes:
            if node not in visited:
                component = []
                kCoreBaseStructuralDiversity.find_component_bfs(k_core_subgraph, node, visited, component)
                if component:
                    components.append(set(component))

        return components

    @staticmethod
    def k_core_decomposition(graph, k):
        """
        K-core decomposition algorithm
        """
        # Handle edge cases
        if not graph or k <= 0:
            if k <= 0 and graph:
                return set(graph.keys())
            return set()

        # Calculate initial degrees for all nodes
        degrees = {node: len(neighbors) for node, neighbors in graph.items()}
        removed = set()
        to_remove = deque()

        # Initialize removal queue - find nodes with degree less than k
        for node, degree in degrees.items():
            if degree < k:
                to_remove.append(node)

        # Iteratively remove nodes with insufficient degree
        while to_remove:
            node = to_remove.popleft()

            # Skip nodes that are already marked as removed
            if node in removed:
                continue

            # Mark current node as removed
            removed.add(node)

            # Update degrees of neighbors
            for neighbor in graph[node]:
                if neighbor not in removed:
                    degrees[neighbor] -= 1
                    # Add neighbor to removal queue if its degree becomes insufficient
                    if degrees[neighbor] < k:
                        to_remove.append(neighbor)

        # Return nodes that are not removed (i.e., k-core nodes)
        return set(graph.keys()) - removed

    @staticmethod
    def find_component_bfs(graph, start, visited, component):
        """Use BFS to find connected component"""
        # Handle edge case: start node already visited
        if start in visited:
            return

        # Handle edge case: start node not in graph
        if start not in graph:
            component.append(start)
            return

        # BFS traversal to find connected component
        queue = deque([start])
        visited.add(start)

        while queue:
            current = queue.popleft()
            component.append(current)

            # Check all neighbors of current node
            if current in graph:
                for neighbor in graph[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)

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
