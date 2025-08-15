#!/usr/bin/env python3
# Auto-generated for 5501746

STUDENT_ID = "5501746"
STUDENT_NAME = "Danlan Zhao"

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
        sd = [0 for _ in range(n)]  # Initialize structural diversity result list

        for v in range(n):
            neighbors = G.adj_list[v]  # Get neighbors of vertex v
            if not neighbors:
                sd[v] = 0
                continue

            # Build a subgraph using neighbors of v (excluding v itself)
            subgraph = SubgraphBuilder.build_neighbor_subgraph(G, neighbors)

            # Compute the number of k-core components in the subgraph
            core_count = CoreCounter.count_k_cores(subgraph, k)
            sd[v] = core_count

        return sd


class SubgraphBuilder:
    @staticmethod
    def build_neighbor_subgraph(G, node_list):
        subgraph = dict()

        for u in node_list:
            subgraph[u] = []

        for u in node_list:
            connections = G.adj_list[u]
            filtered = []
            for v in connections:
                if v in subgraph:
                    filtered.append(v)  # Keep only edges to other neighbors
            subgraph[u].extend(filtered)

        return subgraph


class CoreCounter:
    @staticmethod
    def count_k_cores(adj_dict, threshold):
        #Count how many k-core connected components exist in the given subgraph.
        if not adj_dict:
            return 0

        #Initialize degrees of all nodes
        deg_map = CoreCounter._initialize_degrees(adj_dict)

        #Remove all nodes whose degree is less than k
        removed = CoreCounter._remove_low_degree_nodes(adj_dict, deg_map, threshold)

        #Collect remaining nodes after removals
        remaining_nodes = [u for u in adj_dict if u not in removed]

        if not remaining_nodes:
            return 0

        #Count how many connected components exist in the k-core
        return CoreCounter._count_connected_components(adj_dict, removed)

    @staticmethod
    def _initialize_degrees(adj_dict):
        #Count initial degree of each node in the subgraph.
        deg = {}
        for u in adj_dict:
            deg[u] = len(adj_dict[u])
        return deg

    @staticmethod
    def _remove_low_degree_nodes(adj_dict, deg_map, k):
        to_remove = deque()
        removed = set()

        # Find all nodes with degree less than k
        for u, deg in deg_map.items():
            if deg < k:
                to_remove.append(u)
                removed.add(u)

        # Iteratively remove and update neighbor degrees
        while to_remove:
            u = to_remove.popleft()
            for v in adj_dict[u]:
                if v not in removed:
                    deg_map[v] -= 1
                    if deg_map[v] < k:
                        removed.add(v)
                        to_remove.append(v)

        return removed

    @staticmethod
    def _count_connected_components(adj_dict, removed_set):
        #Use DFS to count the number of connected components in the remaining subgraph.
        visited = set()
        component_count = 0

        for u in adj_dict:
            if u in removed_set or u in visited:
                continue
            component_count += 1
            CoreCounter._dfs(u, adj_dict, visited, removed_set)

        return component_count

    @staticmethod
    def _dfs(start, adj_dict, visited, removed_set):
        #Perform depth-first search on the remaining nodes to mark connected components.
        stack = [start]
        visited.add(start)

        while stack:
            node = stack.pop()
            for neighbor in adj_dict[node]:
                if neighbor not in visited and neighbor not in removed_set:
                    visited.add(neighbor)
                    stack.append(neighbor)
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
