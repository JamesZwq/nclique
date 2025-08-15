#!/usr/bin/env python3
# Auto-generated for 5525102

STUDENT_ID = "5525102"
STUDENT_NAME = "Haoyu Wei"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
               # Get adjacency list representation from the graph
        # Try different possible attributes for graph representation
        if hasattr(G, 'adj_list'):
            adj_list = G.adj_list
        elif hasattr(G, 'adjacency_list'):
            adj_list = G.adjacency_list
        elif hasattr(G, 'graph'):
            adj_list = G.graph
        else:
            # If no standard attribute found, try to construct from edges
            adj_list = [[] for _ in range(n)]
            if hasattr(G, 'edges'):
                for u, v in G.edges:
                    adj_list[u].append(v)
                    adj_list[v].append(u)
            elif hasattr(G, 'edge_list'):
                for u, v in G.edge_list:
                    adj_list[u].append(v)
                    adj_list[v].append(u)
            else:
                # Last resort: assume it has some iterable structure
                try:
                    for u in range(n):
                        for v in G[u]:
                            adj_list[u].append(v)
                except:
                    # If all else fails, return all zeros
                    return sd

        # For each vertex, compute its k-core-based structural diversity
        for v in range(n):
            # Get neighbors of vertex v
            neighbors = set(adj_list[v]) if v < len(adj_list) else set()

            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # Build neighbor-induced subgraph
            neighbor_list = list(neighbors)
            neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbor_list)}

            # Create adjacency list for neighbor-induced subgraph
            adj = [[] for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list):
                for neighbor_v in adj_list[u]:
                    if neighbor_v in neighbor_to_idx and neighbor_v != u:
                        j = neighbor_to_idx[neighbor_v]
                        adj[i].append(j)

            # Find all k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity._find_k_cores(adj, k)
            sd[v] = len(k_cores)
        return sd
    @staticmethod
    def _find_k_cores(adj, k):
        """
        Find all k-cores in a graph represented by adjacency list
        Returns list of k-cores, where each k-core is a set of vertices
        """
        n = len(adj)
        if n == 0:
            return []

        # Compute degrees
        degree = [len(adj[i]) for i in range(n)]

        # Remove vertices with degree < k iteratively
        removed = [False] * n
        queue = deque()

        # Initialize queue with vertices having degree < k
        for i in range(n):
            if degree[i] < k:
                queue.append(i)

        # Remove vertices with degree < k
        while queue:
            u = queue.popleft()
            if removed[u]:
                continue
            removed[u] = True

            # Update degrees of neighbors
            for v in adj[u]:
                if not removed[v]:
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)

        # Find connected components in remaining graph
        remaining = [i for i in range(n) if not removed[i]]
        if not remaining:
            return []

        # Build adjacency list for remaining vertices
        remaining_set = set(remaining)
        remaining_adj = {}
        for u in remaining:
            remaining_adj[u] = []
            for v in adj[u]:
                if v in remaining_set:
                    remaining_adj[u].append(v)

        # Find connected components using DFS
        visited = set()
        components = []

        def dfs(u, component):
            visited.add(u)
            component.add(u)
            for v in remaining_adj[u]:
                if v not in visited:
                    dfs(v, component)

        for u in remaining:
            if u not in visited:
                component = set()
                dfs(u, component)
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
