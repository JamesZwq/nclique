#!/usr/bin/env python3
# Auto-generated for 5514713

STUDENT_ID = "5514713"
STUDENT_NAME = "Hang Wang"

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

        # Compute τ_k(v) for every vertex
        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                sd[v] = 0
                continue

            # Build neighbor-induced subgraph
            neighbor_subgraph = kCoreBaseStructuralDiversity._build_neighbor_subgraph(
                G, v
            )

            # Find all maximal k-cores in neighbor subgraph
            k_cores = kCoreBaseStructuralDiversity._find_maximal_k_cores(
                neighbor_subgraph, k
            )

            # Count number of k-cores
            sd[v] = len(k_cores)

        return sd

    @staticmethod
    def _build_neighbor_subgraph(G, vertex):
        """
        Build the subgraph induced by neighbors
        """
        neighbors = G.adj_list[vertex]
        if not neighbors:
            return {}

        # Map original vertex IDs to local IDs
        vertex_map = {v: i for i, v in enumerate(neighbors)}

        # Build new subgraph
        subgraph = {i: set() for i in range(len(neighbors))}

        for u, i in vertex_map.items():
            for neighbor in G.adj_list[u]:
                if neighbor in vertex_map:
                    j = vertex_map[neighbor]
                    subgraph[i].add(j)

        # Return subgraph with local vertex ID
        return subgraph

    @staticmethod
    def _find_maximal_k_cores(subgraph, k):
        """
        Find all maximal k-cores in the subgraph.
        """
        if not subgraph:
            return []

        # Find all connected components
        components = kCoreBaseStructuralDiversity._find_connected_components(subgraph)

        maximal_k_cores = []

        # Find maximal k-cores
        for component in components:
            # Build subgraph for this component
            component_graph = {v: set() for v in component}
            for v in component:
                for neighbor in subgraph[v]:
                    if neighbor in component:
                        component_graph[v].add(neighbor)

            # Find maximal k-cores in this component
            cores = kCoreBaseStructuralDiversity._find_maximal_k_cores_in_component(
                component_graph, k
            )
            maximal_k_cores.extend(cores)

        return maximal_k_cores

    @staticmethod
    def _find_maximal_k_cores_in_component(graph: dict[int, set[int]], k: int):
        """
        Find all maximal k-cores in a single connected component.
        """
        if not graph:
            return []

        # Queue for vertices to remove
        queue = deque()

        # Initialize queue with vertices having degree < k
        for v in graph:
            if len(graph[v]) < k:
                queue.append(v)

        # Removal process
        while queue:
            v = queue.popleft()
            if v not in graph:
                continue

            # Remove vertex and update its neighbors
            neighbors_to_update = list(graph[v])
            del graph[v]

            for neighbor in neighbors_to_update:
                if neighbor in graph:
                    graph[neighbor].remove(v)
                    # If neighbor now has degree < k, add to queue
                    if len(graph[neighbor]) < k:
                        queue.append(neighbor)

        if not graph:
            return []

        # Each connected component is a maximal k-core now
        return kCoreBaseStructuralDiversity._find_connected_components(graph)

    @staticmethod
    def _find_connected_components(graph):
        """
        Find all connected components in the graph using DFS
        """
        if not graph:
            return []

        visited = set()
        components = []

        for start_vertex in graph:
            # DFS to find connected component
            if start_vertex in visited:
                continue
            component = set()
            stack = [start_vertex]

            while stack:
                v = stack.pop()
                if v not in visited:
                    visited.add(v)
                    component.add(v)

                    # Add unvisited neighbors to stack
                    for neighbor in graph.get(v, []):
                        if neighbor not in visited:
                            stack.append(neighbor)

            if component:
                components.append(component)

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
