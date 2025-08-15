#!/usr/bin/env python3
# Auto-generated for 5474825

STUDENT_ID = "5474825"
STUDENT_NAME = "Shutong Chen"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
import collections
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _get_induced_subgraph(original_graph_adj, vertex_subset_nodes):
        induced_adj = collections.defaultdict(set)
        for u in vertex_subset_nodes:
            for v in original_graph_adj.get(u, set()):
                if v in vertex_subset_nodes:
                    induced_adj[u].add(v)
                    induced_adj[v].add(u)
        return induced_adj

    @staticmethod
    def _find_k_cores_in_subgraph(subgraph_adj, k):
        if not subgraph_adj:
            return []

        current_degrees = collections.defaultdict(int)
        active_vertices = set(subgraph_adj.keys())

        for u in active_vertices:
            current_degrees[u] = len(subgraph_adj[u])

        q = collections.deque([v for v in active_vertices if current_degrees[v] < k])

        while q:
            u = q.popleft()
            if u not in active_vertices:
                continue

            active_vertices.remove(u)

            for v in subgraph_adj[u]:
                if v in active_vertices:
                    current_degrees[v] -= 1
                    if current_degrees[v] < k:
                        q.append(v)

        k_cores = []
        visited_cc = set()

        for v_start in active_vertices:
            if v_start not in visited_cc:
                current_component_nodes = set()
                q_cc = collections.deque([v_start])
                visited_cc.add(v_start)
                current_component_nodes.add(v_start)

                while q_cc:
                    u_cc = q_cc.popleft()
                    for v_cc_neighbor in subgraph_adj.get(u_cc, set()):
                        if v_cc_neighbor in active_vertices and v_cc_neighbor not in visited_cc:
                            visited_cc.add(v_cc_neighbor)
                            current_component_nodes.add(v_cc_neighbor)
                            q_cc.append(v_cc_neighbor)

                is_k_core = True
                if not current_component_nodes:
                    is_k_core = False
                else:
                    for node_in_comp in current_component_nodes:
                        degree_in_comp = 0
                        for neighbor_in_subgraph in subgraph_adj.get(node_in_comp, set()):
                            if neighbor_in_subgraph in current_component_nodes:
                                degree_in_comp += 1
                        if degree_in_comp < k:
                            is_k_core = False
                            break

                if is_k_core:
                    k_cores.append(frozenset(current_component_nodes))

        return k_cores


    @staticmethod
    def process(G, k):
        graph_adj_dict_of_sets = collections.defaultdict(set)
        for u in range(G.vertex_num):
            for v_neighbor in G.adj_list[u]:
                graph_adj_dict_of_sets[u].add(v_neighbor)
                graph_adj_dict_of_sets[v_neighbor].add(u)

        τ = [0] * G.vertex_num

        for v in range(G.vertex_num):
            neighbors_of_v = graph_adj_dict_of_sets.get(v, set())

            if not neighbors_of_v:
                τ[v] = 0
                continue

            nbr_v_subgraph_adj = kCoreBaseStructuralDiversity._get_induced_subgraph(
                graph_adj_dict_of_sets, neighbors_of_v
            )

            k_cores_in_nbr_v = kCoreBaseStructuralDiversity._find_k_cores_in_subgraph(
                nbr_v_subgraph_adj, k
            )

            τ[v] = len(k_cores_in_nbr_v)

        return τ

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
