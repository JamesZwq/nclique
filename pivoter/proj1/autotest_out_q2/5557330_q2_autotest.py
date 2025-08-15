#!/usr/bin/env python3
# Auto-generated for 5557330

STUDENT_ID = "5557330"
STUDENT_NAME = "Xinyuan Li"

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
        structural_diversity = [0] * n

        for vertex in range(n):  # For each vertex, compute its k-core-based structural diversity
            neighbours = G.adj_list[vertex]  # Get all neighbours of the current vertex
            if not neighbours:  # If the vertex has no neighbours, its structural diversity is 0
                structural_diversity[vertex] = 0
                continue

            # Build the subgraph induced by the neighbours
            neighbour_subgraph = {u_node: set() for u_node in neighbours}
            for u_node in neighbours:
                for w_node in G.adj_list[u_node]:
                    if w_node in neighbours:
                        neighbour_subgraph[u_node].add(w_node)

            # Convert sets to lists if needed later
            neighbour_subgraph = {u_node: list(vs) for u_node, vs in neighbour_subgraph.items()}

            # Compute the number of connected components in the k-core of the neighbour-induced subgraph
            structural_diversity[vertex] = kCoreBaseStructuralDiversity._compute_k_core(neighbour_subgraph, k)

        # Optional: comment out for silent mode
        # print(structural_diversity)
        return structural_diversity

    @staticmethod
    def _compute_k_core(neighbour_subgraph, k):
        if not neighbour_subgraph:
            return 0

        # Count degree for each node
        degree_map = {u_node: len(neighbours) for u_node, neighbours in neighbour_subgraph.items()}
        removed_nodes = set()

        # Initialize queue with nodes whose degree < k
        queue_nodes = deque()
        for u_node in neighbour_subgraph:
            if degree_map[u_node] < k:  # Mark nodes with degree less than k for deletion
                queue_nodes.append(u_node)
                removed_nodes.add(u_node)

        # Iteratively remove nodes that violate k-core condition
        while queue_nodes:  # Remove all nodes not satisfying the k-core condition
            u_node = queue_nodes.popleft()
            for vertex in neighbour_subgraph[u_node]:
                if vertex not in removed_nodes:
                    degree_map[vertex] -= 1
                    if degree_map[vertex] < k:
                        removed_nodes.add(vertex)
                        queue_nodes.append(vertex)

        # Remaining nodes form valid k-core(s)
        remaining_nodes = [u_node for u_node in neighbour_subgraph if u_node not in removed_nodes]
        if not remaining_nodes:  # If no nodes remain, there are 0 k-core components
            return 0

        visited_nodes = set()
        k_core_components = 0

        # Count connected components using BFS
        for node in remaining_nodes:  # Count the number of connected components (k-cores)
            if node not in visited_nodes:
                k_core_components += 1
                queue = deque([node])
                visited_nodes.add(node)
                while queue:
                    u_node = queue.popleft()
                    for v_node in neighbour_subgraph[u_node]:
                        if v_node not in removed_nodes and v_node not in visited_nodes:
                            visited_nodes.add(v_node)
                            queue.append(v_node)

        return k_core_components



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
