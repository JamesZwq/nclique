#!/usr/bin/env python3
# Auto-generated for 5498990

STUDENT_ID = "5498990"
STUDENT_NAME = "Yanchuan Xia"

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
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            # Get neighbors of vertex v
            neighbors = G.adj_list[v]

            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # Create neighbor-induced subgraph
            # Map original vertex IDs to local IDs (0, 1, 2, ...)
            vertex_map = {neighbors[i]: i for i in range(len(neighbors))}
            local_adj = [[] for _ in range(len(neighbors))]

            # Build adjacency list for neighbor-induced subgraph
            for i, u in enumerate(neighbors):
                for w in G.adj_list[u]:
                    if w in vertex_map:  # w is also a neighbor of v
                        local_adj[i].append(vertex_map[w])

            # Find all k-cores in the neighbor-induced subgraph
            sd[v] = kCoreBaseStructuralDiversity._count_k_cores(local_adj, k)

        return sd

    @staticmethod
    def _count_k_cores(adj_list, k):
        """
        Count the number of k-cores in the given graph
        """
        n = len(adj_list)
        if n == 0:
            return 0

        # Find all k-cores using iterative approach
        visited_global = [False] * n
        k_core_count = 0

        for start in range(n):
            if visited_global[start]:
                continue

            # Get the connected component
            component = []
            queue = deque([start])
            visited_comp = [False] * n
            visited_comp[start] = True
            visited_global[start] = True

            while queue:
                v = queue.popleft()
                component.append(v)

                for u in adj_list[v]:
                    if not visited_comp[u]:
                        visited_comp[u] = True
                        visited_global[u] = True
                        queue.append(u)

            # Find all k-cores within this component
            k_cores = kCoreBaseStructuralDiversity._find_all_k_cores_in_component(adj_list, component, k)
            k_core_count += len(k_cores)

        return k_core_count



    @staticmethod
    def _find_all_k_cores_in_component(adj_list, component, k):
        """
        Find all maximal k-cores within a given component
        """
        if len(component) == 0:
            return []

        # Create subgraph for this component
        comp_set = set(component)

        # Calculate initial degrees within the component
        degrees = {}
        for v in component:
            degree = 0
            for u in adj_list[v]:
                if u in comp_set:
                    degree += 1
            degrees[v] = degree

        # Remove vertices with degree < k iteratively to get the k-core
        queue = deque()
        removed = set()

        # Initialize queue with vertices having degree < k
        for v in component:
            if degrees[v] < k:
                queue.append(v)
                removed.add(v)

        # Iteratively remove vertices
        while queue:
            v = queue.popleft()

            for u in adj_list[v]:
                if u in comp_set and u not in removed:
                    degrees[u] -= 1
                    if degrees[u] < k and u not in removed:
                        queue.append(u)
                        removed.add(u)

        # Get the remaining vertices (they form the k-core subgraph)
        k_core_vertices = [v for v in component if v not in removed]

        if len(k_core_vertices) == 0:
            return []

        # Now find all maximal k-cores within this k-core subgraph
        # A maximal k-core is a connected component within the k-core subgraph
        k_cores = []
        visited = set()

        for v in k_core_vertices:
            if v in visited:
                continue

            # Find connected component of k-core vertices starting from v
            core_component = []
            queue = deque([v])
            visited.add(v)

            while queue:
                curr = queue.popleft()
                core_component.append(curr)

                for neighbor in adj_list[curr]:
                    if neighbor in k_core_vertices and neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)

            if len(core_component) > 0:
                k_cores.append(core_component)

        return k_cores

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
