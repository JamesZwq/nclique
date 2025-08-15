#!/usr/bin/env python3
# Auto-generated for 5448528

STUDENT_ID = "5448528"
STUDENT_NAME = "Jingting Zhou"

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

        # For each vertex, compute its k-core
        for v in range(n):
            # Get neighbors of v
            neighbors = set(G.adj_list[v])

            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # Create neighbor-induced subgraph
            neighbor_induced = kCoreBaseStructuralDiversity._create_neighbor_induced_subgraph(G, v, neighbors)

            # Find k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity._find_k_cores(neighbor_induced, k)

            # Count the number of k-cores
            sd[v] = len(k_cores)

        return sd

    @staticmethod
    def _create_neighbor_induced_subgraph(G, vertex, neighbors):
        """
        Create the neighbor-induced subgraph G[N(v)] for vertex v.
        Returns adjacency list representation of the subgraph.
        """
        neighbor_list = list(neighbors)
        vertex_to_index = {v: i for i, v in enumerate(neighbor_list)}

        # Build adjacency list for neighbor-induced subgraph
        subgraph_adj = [[] for _ in range(len(neighbor_list))]

        for i, u in enumerate(neighbor_list):
            for neighbor in G.adj_list[u]:
                if neighbor in neighbors and neighbor != u:
                    j = vertex_to_index[neighbor]
                    subgraph_adj[i].append(j)

        # Remove duplicates and sort
        for i in range(len(subgraph_adj)):
            subgraph_adj[i] = sorted(list(set(subgraph_adj[i])))

        return subgraph_adj, neighbor_list

    @staticmethod
    def _find_k_cores(subgraph_data, k):

        adj_list, original_vertices = subgraph_data
        n = len(adj_list)

        if n == 0:
            return []

        # Compute k-core decomposition

        remaining = set(range(n))
        degrees = [len(adj_list[i]) for i in range(n)]

        changed = True
        while changed:
            changed = False
            to_remove = []

            for v in remaining:
                if degrees[v] < k:
                    to_remove.append(v)

            for v in to_remove:
                remaining.remove(v)
                changed = True

                # Update degrees of neighbors
                for neighbor in adj_list[v]:
                    if neighbor in remaining:
                        degrees[neighbor] -= 1

        # Find connected components in remaining vertices
        if not remaining:
            return []

        # Build adjacency list for remaining vertices
        remaining_adj = {v: [] for v in remaining}
        for v in remaining:
            for neighbor in adj_list[v]:
                if neighbor in remaining:
                    remaining_adj[v].append(neighbor)

        # Find connected components using DFS
        visited = set()
        k_cores = []

        for start_vertex in remaining:
            if start_vertex not in visited:
                # DFS
                component = []
                stack = [start_vertex]

                while stack:
                    v = stack.pop()
                    if v not in visited:
                        visited.add(v)
                        component.append(original_vertices[v])

                        for neighbor in remaining_adj[v]:
                            if neighbor not in visited:
                                stack.append(neighbor)

                if component:
                    k_cores.append(sorted(component))

        return k_cores

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
