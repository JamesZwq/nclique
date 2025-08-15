#!/usr/bin/env python3
# Auto-generated for 5609829

STUDENT_ID = "5609829"
STUDENT_NAME = "Song Huang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################



import sys



class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _create_neighbor_subgraph(neighbors_list, full_adj_list):

        neighbor_set = set(neighbors_list)
        subgraph_adj = {node: [] for node in neighbors_list}
        subgraph_degrees = {node: 0 for node in neighbors_list}

        for u_node in neighbors_list:
            for v_node in full_adj_list[u_node]:
                if v_node in neighbor_set:
                    if u_node < v_node:
                        subgraph_adj[u_node].append(v_node)
                        subgraph_adj[v_node].append(u_node)
                        subgraph_degrees[u_node] += 1
                        subgraph_degrees[v_node] += 1

        return subgraph_adj, subgraph_degrees

    @staticmethod
    def _get_core_nodes(subgraph_adj, subgraph_degrees, k_val):

        degrees = subgraph_degrees.copy()
        nodes_in_core = set(degrees.keys())


        removal_queue = deque([
            node for node, deg in degrees.items() if deg < k_val
        ])

        while removal_queue:
            node = removal_queue.popleft()
            if node not in nodes_in_core:
                continue

            nodes_in_core.remove(node)


            for neighbor in subgraph_adj[node]:
                if neighbor in nodes_in_core:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] == k_val - 1:
                        removal_queue.append(neighbor)

        return nodes_in_core

    @staticmethod
    def _dfs_recursive_traversal(node, adj, core_nodes, visited):

        visited.add(node)
        for neighbor in adj[node]:
            if neighbor in core_nodes and neighbor not in visited:
                kCoreBaseStructuralDiversity._dfs_recursive_traversal(
                    neighbor, adj, core_nodes, visited
                )

    @staticmethod
    def process(G, k):

        final_tau_values = [0] * G.vertex_num

        for v_idx in range(G.vertex_num):
            neighbors = G.adj_list[v_idx]


            if len(neighbors) <= k:
                continue


            subgraph_adj, subgraph_degrees = kCoreBaseStructuralDiversity._create_neighbor_subgraph(
                neighbors, G.adj_list
            )


            core_nodes = kCoreBaseStructuralDiversity._get_core_nodes(
                subgraph_adj, subgraph_degrees, k
            )


            if not core_nodes:
                num_components = 0
            else:
                visited_nodes = set()
                num_components = 0
                for node in core_nodes:
                    if node not in visited_nodes:
                        num_components += 1

                        kCoreBaseStructuralDiversity._dfs_recursive_traversal(
                            node, subgraph_adj, core_nodes, visited_nodes
                        )

            final_tau_values[v_idx] = num_components

        return final_tau_values

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
