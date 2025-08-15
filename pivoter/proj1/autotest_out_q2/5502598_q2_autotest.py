#!/usr/bin/env python3
# Auto-generated for 5502598

STUDENT_ID = "5502598"
STUDENT_NAME = "Yufei Liu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
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

        # For each vertex, compute k-core based structural diversity
        for v in range(n):
            neighbors = set(G.adj_list[v])
            if len(neighbors) == 0:
                tau[v] = 0
                continue

            # Build neighbor-induced subgraph
            neighbor_graph = kCoreBaseStructuralDiversity._build_neighbor_subgraph(G, neighbors)

            # Find k-cores in the neighbor-induced subgraph
            k_core_components = kCoreBaseStructuralDiversity._find_k_cores(neighbor_graph, k)

            tau[v] = len(k_core_components)

        return tau

    @staticmethod
    def _build_neighbor_subgraph(G, neighbors):
        """
        Build induced subgraph from neighbor vertices
        Returns adjacency list representation
        """
        neighbor_list = list(neighbors)
        neighbor_to_idx = {v: i for i, v in enumerate(neighbor_list)}

        # Build adjacency list for neighbor subgraph
        subgraph = [[] for _ in range(len(neighbor_list))]

        for i, u in enumerate(neighbor_list):
            for v in G.adj_list[u]:
                if v in neighbors:
                    j = neighbor_to_idx[v]
                    subgraph[i].append(j)

        return subgraph, neighbor_list

    @staticmethod
    def _find_k_cores(neighbor_graph_info, k):
        """
        Find all k-cores (connected components where all vertices have degree >= k)
        """
        subgraph, original_vertices = neighbor_graph_info
        n = len(subgraph)

        if n == 0:
            return []

        # Step 1: Remove vertices with degree < k iteratively
        active = [True] * n
        degrees = [len(adj) for adj in subgraph]
        queue = deque()

        # Initialize queue with vertices having degree < k
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)
                active[i] = False

        # Iteratively remove vertices with degree < k
        while queue:
            u = queue.popleft()
            for v in subgraph[u]:
                if active[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)
                        active[v] = False

        # Step 2: Find connected components in remaining vertices
        remaining_vertices = [i for i in range(n) if active[i]]

        if not remaining_vertices:
            return []

        # Build adjacency list for remaining vertices only
        remaining_adj = defaultdict(list)
        for u in remaining_vertices:
            for v in subgraph[u]:
                if active[v]:
                    remaining_adj[u].append(v)

        # Find connected components using DFS
        visited = [False] * n
        components = []

        for start in remaining_vertices:
            if not visited[start]:
                component = []
                stack = [start]

                while stack:
                    u = stack.pop()
                    if not visited[u]:
                        visited[u] = True
                        component.append(original_vertices[u])

                        for v in remaining_adj[u]:
                            if not visited[v]:
                                stack.append(v)

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
