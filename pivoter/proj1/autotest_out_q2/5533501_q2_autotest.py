#!/usr/bin/env python3
# Auto-generated for 5533501

STUDENT_ID = "5533501"
STUDENT_NAME = "Yinjie Ma"

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

        def _count_connected_components(adj_dict):
            """
            Compute connected components in a graph
            """
            if not adj_dict:
                return 0

            visited_nodes = set()
            component_count = 0

            for vertex in adj_dict:
                if vertex not in visited_nodes:
                    # Use breadth-first search (BFS) to traverse connected components
                    queue = deque([vertex])
                    visited_nodes.add(vertex)

                    while queue:
                        current_node = queue.popleft()
                        for neighbor in adj_dict[current_node]:
                            if neighbor not in visited_nodes:
                                visited_nodes.add(neighbor)
                                queue.append(neighbor)

                    component_count += 1

            return component_count

        def _compute_k_cores(adj_dict, core_value):
            """
            Count the number of k-cores in the graph
            """
            if not adj_dict or core_value < 0:
                return 0

            if core_value == 0:
                return _count_connected_components(adj_dict)

            graph_copy = {vertex: set(neighbors) for vertex, neighbors in adj_dict.items()}

            updated = True
            while updated:
                updated = False
                nodes_to_remove = []

                # Find nodes with degree less than k
                for vertex in graph_copy:
                    if len(graph_copy[vertex]) < core_value:
                        nodes_to_remove.append(vertex)

                # Delete nodes with degree less than k
                for vertex in nodes_to_remove:
                    for neighbor in graph_copy[vertex]:
                        graph_copy[neighbor].discard(vertex)
                    del graph_copy[vertex]
                    updated = True

            if not graph_copy:
                return 0

            visited_nodes = set()
            component_count = 0

            # Calculate the number of connected components of the remaining subgraph
            for vertex in graph_copy:
                if vertex not in visited_nodes:
                    queue = deque([vertex])
                    visited_nodes.add(vertex)

                    while queue:
                        current_node = queue.popleft()
                        for neighbor in graph_copy[current_node]:
                            if neighbor not in visited_nodes:
                                visited_nodes.add(neighbor)
                                queue.append(neighbor)

                    component_count += 1

            return component_count

        for vertex in range(n):
            # Get the neighbors of a vertex
            neighbors = G.adj_list[vertex]
            if len(neighbors) == 0:
                sd[vertex] = 0
                continue

            # Constructing neighbor-induced subgraphs
            neighbor_set = set(neighbors)

            # Create adjacency list of neighbor-induced subgraph
            subgraph_adj = {}
            for neighbor in neighbors:
                subgraph_adj[neighbor] = []
                for adjacent in G.adj_list[neighbor]:
                    if adjacent in neighbor_set:
                        subgraph_adj[neighbor].append(adjacent)

            # Calculate k-cores in a subgraph
            sd[vertex] = _compute_k_cores(subgraph_adj, k)

        return sd

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
