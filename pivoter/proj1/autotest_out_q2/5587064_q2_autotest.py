#!/usr/bin/env python3
# Auto-generated for 5587064

STUDENT_ID = "5587064"
STUDENT_NAME = "Zisen Wang"

# ======= 学生代码 =======
from collections import deque
import sys


sys.setrecursionlimit(2000000)

################################################################################
# You can import any Python Standard Library modules~
################################################################################
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _find_k_core_nodes(subgraph_adj, k):

        node_count = len(subgraph_adj)
        if node_count == 0:
            return set()

        degrees = {i: len(neighbors) for i, neighbors in enumerate(subgraph_adj)}

       # Queue for nodes to be removed
        removal_queue = deque([i for i, deg in degrees.items() if deg < k])

        removed_nodes = set(removal_queue)

        while removal_queue:
            node_to_remove = removal_queue.popleft()

            for neighbor in subgraph_adj[node_to_remove]:
                if neighbor not in removed_nodes:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        removed_nodes.add(neighbor)
                        removal_queue.append(neighbor)

        # The core nodes are those that were never removed.
        all_nodes = set(range(node_count))
        core_nodes = all_nodes - removed_nodes
        return core_nodes

    @staticmethod
    def _dfs_component_search(node, subgraph_adj, core_nodes, visited):

        visited.add(node)
        for neighbor in subgraph_adj[node]:
            if neighbor in core_nodes and neighbor not in visited:
                kCoreBaseStructuralDiversity._dfs_component_search(neighbor, subgraph_adj, core_nodes, visited)

    @staticmethod
    def _count_components(subgraph_adj, core_nodes):

        if not core_nodes:
            return 0

        component_count = 0
        visited = set()

        for node in core_nodes:
            if node not in visited:
                component_count += 1
                kCoreBaseStructuralDiversity._dfs_component_search(node, subgraph_adj, core_nodes, visited)

        return component_count

    @staticmethod
    def process(G, k):

        num_vertices = G.vertex_num
        tau_vector = [0] * num_vertices

        for center_node in range(num_vertices):

            neighbors = G.adj_list[center_node]
            if not neighbors:
                continue


            local_node_map = {neighbor_id: i for i, neighbor_id in enumerate(neighbors)}


            subgraph_adj = [[] for _ in range(len(neighbors))]
            for original_id, local_idx in local_node_map.items():
                for neighbor_of_neighbor in G.adj_list[original_id]:
                    if neighbor_of_neighbor in local_node_map:
                        subgraph_adj[local_idx].append(local_node_map[neighbor_of_neighbor])


            core_nodes = kCoreBaseStructuralDiversity._find_k_core_nodes(subgraph_adj, k)


            num_k_cores = kCoreBaseStructuralDiversity._count_components(subgraph_adj, core_nodes)

            tau_vector[center_node] = num_k_cores

        return tau_vector

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
