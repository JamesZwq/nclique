#!/usr/bin/env python3
# Auto-generated for 5098800

STUDENT_ID = "5098800"
STUDENT_NAME = "Shengyun Yin"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _get_induced_subgraph(graph, vertices_subset):
        """
        Creates an induced subgraph from a given graph and a subset of its vertices.
        The new graph will have re-indexed vertices starting from 0.
        """
        if not vertices_subset:
            return UndirectedUnweightedGraph([(0, 0)])

        # Create mapping from original ID to new ID
        original_to_new_id = {v_id: i for i, v_id in enumerate(list(vertices_subset))}

        subgraph_vertex_num = len(vertices_subset)
        subgraph_edge_list_for_constructor = []
        subgraph_edge_count = 0

        # Build edges for the new UndirectedUnweightedGraph constructor
        # Iterate over original vertices in the subset
        for u_orig in list(vertices_subset):
            # Iterate over neighbors of u_orig in the original graph
            for v_orig in graph.adj_list[u_orig]:
                # If v_orig is also in the subset, add edges
                if v_orig in vertices_subset and u_orig < v_orig:
                    u_new = original_to_new_id[u_orig]
                    v_new = original_to_new_id[v_orig]
                    subgraph_edge_list_for_constructor.append((u_new, v_new))
                    subgraph_edge_count += 1

        constructor_input = [(subgraph_vertex_num, subgraph_edge_count)] + subgraph_edge_list_for_constructor
        return UndirectedUnweightedGraph(constructor_input)

    @staticmethod
    def _get_k_core_and_components(graph, k):
        """
        Computes the k-core of the given graph and returns the number of connected components within it.
        """
        n_sub = graph.vertex_num
        if n_sub == 0:
            return 0 # No vertices, no k-core, 0 components

        # Step 1: Compute the k-core using iterative pruning
        current_degrees = [len(graph.adj_list[i]) for i in range(n_sub)]

        # Maintain a set of vertices currently in the core candidate
        core_vertices_set = set(range(n_sub))

        # Queue for vertices to be removed
        q = deque()

        # Initialize queue with vertices having degree < k
        for i in range(n_sub):
            if current_degrees[i] < k:
                q.append(i)
                if i in core_vertices_set: # Ensure it's not already removed from set
                    core_vertices_set.remove(i)

        # Iterative pruning
        while q:
            u = q.popleft()

            # For each neighbor of u in the ORIGINAL graph (which is graph here)
            for v in graph.adj_list[u]:
                if v in core_vertices_set: # If neighbor is still a candidate for the core
                    current_degrees[v] -= 1 # Decrement its degree
                    if current_degrees[v] < k:
                        q.append(v)
                        core_vertices_set.remove(v) # Remove from core candidates

        # Step 2: Count connected components in the k-core
        num_components = 0
        visited = set()

        for v_core in core_vertices_set:
            if v_core not in visited:
                num_components += 1
                q_bfs = deque([v_core])
                visited.add(v_core)

                while q_bfs:
                    curr_node = q_bfs.popleft()
                    # Only traverse to neighbors that are also in the k-core and haven't been visited
                    for neighbor in graph.adj_list[curr_node]:
                        if neighbor in core_vertices_set and neighbor not in visited:
                            visited.add(neighbor)
                            q_bfs.append(neighbor)

        return num_components

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
            neighbors = G.adj_list[v]

            if not neighbors:
                sd[v] = 0
                continue

            # Create the neighbour-induced subgraph for vertex v
            nbr_v_graph = kCoreBaseStructuralDiversity._get_induced_subgraph(G, set(neighbors))

            # Calculate the number of k-cores in the neighbor-induced subgraph
            sd[v] = kCoreBaseStructuralDiversity._get_k_core_and_components(nbr_v_graph, k)

        return sd

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
