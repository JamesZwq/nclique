#!/usr/bin/env python3
# Auto-generated for 5520110

STUDENT_ID = "5520110"
STUDENT_NAME = "Yiling Jiang"

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
        Computes τ_k(v) for each vertex v in the graph G.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The input graph in adjacency list format.
        k : int
            The k-core threshold.

        Returns
        -------
        List[int]
            A list where the i-th element is τ_k(i),
            i.e., the number of k-cores in the neighbour-induced subgraph of vertex i.
        """
        n = G.vertex_num
        sd = [0] * n  # Output array storing τ_k(v) for each vertex

        for v in range(n):
            neighbors = G.adj_list[v]
            nbr_set = set(neighbors)

            # Step 1: Build the neighbour-induced subgraph of v
            # Only include edges between neighbors of v
            subgraph = {u: [] for u in nbr_set}
            for u in nbr_set:
                for w in G.adj_list[u]:
                    if w in nbr_set:
                        subgraph[u].append(w)

            # Step 2: Prune the subgraph to its k-core
            core_nodes = kCoreBaseStructuralDiversity.k_core_prune(subgraph, k)

            # Step 3: Count connected components in the pruned subgraph
            comp_count = kCoreBaseStructuralDiversity.connected_components(
                subgraph, core_nodes, G.vertex_num
            )

            sd[v] = comp_count  # τ_k(v)

        return sd

    @staticmethod
    def k_core_prune(subgraph, k):
        """
        Removes nodes from the subgraph that do not satisfy the k-core condition.

        Parameters
        ----------
        subgraph : dict[int, List[int]]
            The neighbour-induced subgraph represented as an adjacency list.
        k : int
            The minimum degree requirement for nodes in the k-core.

        Returns
        -------
        Set[int]
            The set of nodes remaining in the k-core of the subgraph.
        """
        # Initialize degrees
        degree = {u: len(subgraph[u]) for u in subgraph}
        # Queue for nodes with degree less than k
        queue = deque([u for u in subgraph if degree[u] < k])
        removed = set()

        # Iteratively remove nodes with degree < k and update their neighbors
        while queue:
            u = queue.popleft()
            removed.add(u)
            for v in subgraph[u]:
                if v not in removed:
                    degree[v] -= 1
                    if degree[v] == k - 1:
                        queue.append(v)

        return set(subgraph.keys()) - removed  # Nodes that survived pruning

    @staticmethod
    def connected_components(subgraph, nodes, total_nodes):
        """
        Counts the number of connected components in a subgraph.

        Parameters
        ----------
        subgraph : dict[int, List[int]]
            The adjacency list of the subgraph.
        nodes : Set[int]
            The nodes to include in the component search (after k-core pruning).
        total_nodes : int
            Total number of nodes in the original graph (for visited array sizing).

        Returns
        -------
        int
            The number of connected components among the given nodes.
        """
        visited = [False] * total_nodes
        count = 0

        for u in nodes:
            if visited[u]:
                continue
            count += 1
            queue = deque([u])
            visited[u] = True

            # Perform BFS to mark all reachable nodes in this component
            while queue:
                s = queue.popleft()
                for v in subgraph[s]:
                    if v in nodes and not visited[v]:
                        visited[v] = True
                        queue.append(v)

        return count

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
