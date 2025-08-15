#!/usr/bin/env python3
# Auto-generated for 5607223

STUDENT_ID = "5607223"
STUDENT_NAME = "Ruihan Liu"

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

        def find_k_core_v(graph_adj_map, k_val):
            if not graph_adj_map:
                return set()

            current_degrees = {v: len(neighbors) for v, neighbors in graph_adj_map.items()}
            q = deque([v for v, deg in current_degrees.items() if deg < k_val])
            a_vertices = set(graph_adj_map.keys()) #dynamic track of points in the graph

            while q:
                u = q.popleft()
                if u not in a_vertices:
                    continue # already deleted

                a_vertices.remove(u)
                for v in graph_adj_map[u]:
                    if v in a_vertices:
                        current_degrees[v] -= 1
                        if current_degrees[v] < k_val: # smaller than k deg
                            q.append(v)

            return a_vertices

        def count_cc(graph_adj_map, component_vertices_subset=None):
            if not graph_adj_map:
                return 0
            if component_vertices_subset is None:
                nodes_to_visit = set(graph_adj_map.keys())
            else:
                nodes_to_visit = set(component_vertices_subset)
            if not nodes_to_visit:
                return 0
            visited = set() # for visited nodes
            num_components = 0 # count connect component
            for start_node in sorted(list(nodes_to_visit)):
                if start_node not in visited and start_node in nodes_to_visit:
                    num_components += 1 # new cc
                    q = deque([start_node])
                    visited.add(start_node)
                    # BFS through all nodes
                    while q:
                        u = q.popleft()
                        for v in graph_adj_map.get(u, []):
                            if v in nodes_to_visit and v not in visited:
                                visited.add(v)
                                q.append(v)
            return num_components

        # start processing
        for v_id in range(n):
            neighbors = G.adj_list[v_id]

            if not neighbors: # no neighbor, skip
                sd[v_id] = 0
                continue

            nbr_v_adj = {}
            for u in neighbors: # every neighbor is a key
                nbr_v_adj[u] = []

            for u in neighbors:
                for potential_neighbor in G.adj_list[u]:
                    if potential_neighbor in nbr_v_adj:
                        nbr_v_adj[u].append(potential_neighbor)

            for u in nbr_v_adj:
                nbr_v_adj[u] = list(set(nbr_v_adj[u]))

            k_core_v = find_k_core_v(nbr_v_adj, k)
            num_components = count_cc(nbr_v_adj, k_core_v)
            sd[v_id] = num_components

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
