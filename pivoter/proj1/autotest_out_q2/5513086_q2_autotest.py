#!/usr/bin/env python3
# Auto-generated for 5513086

STUDENT_ID = "5513086"
STUDENT_NAME = "Linhui Liu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    @staticmethod
    def process(G, k):
        """
        Compute the k-core-based structural diversity τ_k(v) for every node v in the graph.
        τ_k(v) is defined as the number of connected components in the k-core of the
        neighbor-induced subgraph of v.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The input undirected, unweighted graph.
        k : int
            Minimum degree threshold for k-core.

        Returns
        -------
        List[int]
            τ_k(v) values for all nodes in the graph.
        """
        n = G.vertex_num
        adj = G.adj_list

        # Step 1: Global core number preprocessing
        # Use Batagelj's algorithm to compute the core number of each node.
        # This allows us to skip nodes and neighbors that cannot be in any k-core.
        core = kCoreBaseStructuralDiversity.compute_core_number(G)

        tau_values = [0] * n  # Initialize τ_k(v) for all nodes

        for v in range(n):
            # Pruning 1: Skip nodes with degree < k or core number < k
            if len(adj[v]) < k or core[v] < k:
                tau_values[v] = 0
                continue

            # Step 2: Neighbor filtering
            # Only keep neighbors with core number >= k.
            neighbors = [u for u in adj[v] if core[u] >= k]
            if len(neighbors) < k:
                # If fewer than k neighbors remain, no k-core is possible.
                tau_values[v] = 0
                continue
            neighbor_set = set(neighbors)

            # Step 3: Build neighbor-induced subgraph using a hash set for fast intersection
            # For each neighbor, only keep edges to other neighbors.
            sub_adj = {u: [] for u in neighbors}
            for u in neighbors:
                for w in adj[u]:
                    if w in neighbor_set:
                        sub_adj[u].append(w)

            # Step 4: Extract k-core from the neighbor-induced subgraph
            # Use queue-based peeling to iteratively remove nodes with degree < k.
            k_core_adj = kCoreBaseStructuralDiversity.compute_k_core(sub_adj, k)
            if not k_core_adj:
                # If the k-core is empty, τ_k(v) = 0.
                tau_values[v] = 0
                continue

            # Step 5: Count the number of connected components in the k-core
            tau_values[v] = kCoreBaseStructuralDiversity.count_components(k_core_adj)

        return tau_values

    @staticmethod
    def compute_core_number(G):
        """
        Compute the core number of each node using Batagelj's algorithm.
        Complexity: O(n + m)

        Parameters
        ----------
        G : UndirectedUnweightedGraph

        Returns
        -------
        List[int]
            Core numbers for all nodes.
        """
        from collections import defaultdict
        n = G.vertex_num
        deg = [len(G.adj_list[v]) for v in range(n)]
        core = [0] * n
        bucket = defaultdict(list)

        # Initialize degree buckets
        for v in range(n):
            bucket[deg[v]].append(v)

        curr_deg = 0
        processed = 0
        while processed < n:
            # Find the smallest non-empty bucket
            while not bucket[curr_deg]:
                curr_deg += 1
            v = bucket[curr_deg].pop()
            core[v] = curr_deg
            processed += 1

            # Decrease the degree of neighbors
            for u in G.adj_list[v]:
                if deg[u] > curr_deg:
                    bucket[deg[u]].remove(u)
                    deg[u] -= 1
                    bucket[deg[u]].append(u)
        return core

    @staticmethod
    def compute_k_core(adj_list, k):
        """
        Extract the k-core of a subgraph using queue-based peeling.

        Parameters
        ----------
        adj_list : dict[int, list[int]]
            Adjacency list of the subgraph.
        k : int
            Minimum degree threshold.

        Returns
        -------
        dict[int, list[int]]
            Pruned adjacency list representing the k-core.
        """
        degrees = {u: len(neigh) for u, neigh in adj_list.items()}
        queue = deque([u for u in adj_list if degrees[u] < k])

        # Iteratively remove nodes with degree < k
        while queue:
            u = queue.popleft()
            for v in adj_list[u]:
                if v in degrees:
                    degrees[v] -= 1
                    if degrees[v] == k - 1:
                        queue.append(v)
            del degrees[u]
            adj_list[u] = []

        # Return the pruned subgraph containing only nodes with degree >= k
        return {u: [w for w in neigh if w in degrees] for u, neigh in adj_list.items() if u in degrees}

    @staticmethod
    def count_components(adj_list):
        """
        Count the number of connected components in a graph using DFS.

        Parameters
        ----------
        adj_list : dict[int, list[int]]
            Adjacency list of the k-core subgraph.

        Returns
        -------
        int
            Number of connected components.
        """
        visited = set()
        comp_count = 0

        for node in adj_list:
            if node not in visited:
                comp_count += 1
                stack = [node]
                visited.add(node)
                # Explicit stack DFS to avoid recursion limits
                while stack:
                    cur = stack.pop()
                    for w in adj_list[cur]:
                        if w not in visited:
                            visited.add(w)
                            stack.append(w)
        return comp_count

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
