#!/usr/bin/env python3
# Auto-generated for 5506082

STUDENT_ID = "5506082"
STUDENT_NAME = "Qian Lu"

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

        def get_neighbors(G, v):
            """
            Returns the neighbors of a vertex v in the graph.

            Parameters
            ----------
            G : UndirectedUnweightedGraph
            v : int

            Returns
            -------
            List[int]
            """
            return G.adj_list[v]

        def compute_global_core_number(G):
            """
            Computes the core number for each vertex in the graph using

            Parameters
            ----------
            G : UndirectedUnweightedGraph

            Returns
            -------
            core[v] = max core number of vertex v
            """
            # Initial the degree of each vertex
            degree = [len(get_neighbors(G, v)) for v in range(G.vertex_num)]
            core = degree[:]
            order = []

            # Create a bucket of degrees
            # bin[d]: Vertex bucket with degree d
            bin = [[] for _ in range(max(degree)+1)]
            for v, d in enumerate(degree):
                bin[d].append(v)

            # core decomposition in ascending order
            for d in range(len(bin)):
                while bin[d]:
                    v = bin[d].pop()
                    core[v] = d
                    order.append(v)
                    for u in get_neighbors(G, v):
                        if core[u] > d:
                            bin[core[u]].remove(u)
                            core[u] -= 1
                            bin[core[u]].append(u)
            return core
        
        def build_subgraph(G, valid_nodes):
            """
            Constructs a subgraph from G containing only valid nodes
            and edges between them.

            Parameters
            ----------
            G : UndirectedUnweightedGraph
            valid_nodes : Set[int]

            Returns
            -------
            Dict[int, Set[int]]  # adjacency list of the subgraph
            """
            subG = defaultdict(set)
            for u in valid_nodes:
                for v in get_neighbors(G, u):
                    if v in valid_nodes:
                        subG[u].add(v)
            return subG

        def get_connected_components(subG):
            """
            Finds all connected components in an undirected subgraph.

            Parameters
            ----------
            subG : Dict[int, Set[int]]

            Returns
            -------
            List[List[int]]  # list of components, each as a list of vertices
            """
            visited = set()
            components = []
            for vertex in subG:
                if vertex not in visited:
                    component = []
                    queue = deque([vertex])
                    while queue:
                        u = queue.popleft()
                        if u in visited:
                            continue
                        visited.add(u)
                        component.append(u)
                        for v in subG[u]:
                            if v not in visited:
                                queue.append(v)
                    components.append(component)
            return components

        def prune_k_core(subG, k):
            """
            Performs k-core pruning on a subgraph by iteratively removing
            nodes with degree < k.

            Parameters
            ----------
            subG : Dict[int, Set[int]]
            k : int

            Returns
            -------
            Dict[int, Set[int]]  # pruned subgraph where all nodes have degree ≥ k
            """
            degree = {u: len(subG[u]) for u in subG}
            queue = deque([u for u in subG if degree[u] < k])
            while queue:
                u = queue.popleft()
                # Update the degrees of u's neighbor
                for v in subG[u]:
                    if v in subG and u in subG[v]:
                        subG[v].remove(u)
                        if len(subG[v]) == k - 1:
                            queue.append(v)
                del subG[u]
            return subG

        def count_components(subG):
            """
            Counts the number of connected components in the given subgraph.

            Parameters
            ----------
            subG : Dict[int, Set[int]]

            Returns
            -------
            int  # number of connected components
            """
            visited = set()
            count = 0
            for u in subG:
                if u not in visited:
                    count += 1
                    queue = deque([u])
                    while queue:
                        v = queue.popleft()
                        if v in visited:
                            continue
                        visited.add(v)
                        for nbr in subG[v]:
                            if nbr not in visited:
                                queue.append(nbr)
            return count

        # ---- Main procedure ----
        n = G.vertex_num
        τ = [0] * n
        core = compute_global_core_number(G)

        # All vertices with cores less than k+1 are set to 0 in τ
        valid_nodes = {v for v in range(n) if core[v] >= k + 1}
        if not valid_nodes:
            return τ

        # Build a subgraph with k+1 core vertices
        kplus_subgraph = build_subgraph(G, valid_nodes)
        components = get_connected_components(kplus_subgraph)

        for component in components:
            unique_component = set(component)
            for v in component:
                neighbor_set = {u for u in get_neighbors(G, v) if u in unique_component and core[u] >= k}
                if not neighbor_set:
                    τ[v] = 0
                    continue

                # Build neighbor-induced subgraph
                induced_subgraph = defaultdict(set)
                for u in neighbor_set:
                    for w in get_neighbors(G, u):
                        if w in neighbor_set:
                            induced_subgraph[u].add(w)

                # Prune the subgraph to keep only k-core
                pruned_subgraph = prune_k_core(induced_subgraph, k)
                τ[v] = count_components(pruned_subgraph)

        return τ

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
