#!/usr/bin/env python3
# Auto-generated for 5505423

STUDENT_ID = "5505423"
STUDENT_NAME = "Zhonghao Tong"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Core Number Calculation
             core_num = compute_core_numbers(G)
             Time Complexity: O(n + m)
             Uses a flat array for k-core decomposition.

        Construct (k+1)-Core Adjacency List
            adj_k1core = coreAdj(G, core_num, k + 1)
            Time Complexity: O(n + m)
            Traverses all nodes and edges.

        Build Subgraphs and Compute k-Core Counts for Each Node
        Per iteration:
            buildSubgraph(...) → O(n′ + m′)
            kCoreCount(...) → O(n′ + m′)

        Worst-case: Up to n iterations (one per node)
        Total Time Complexity (Worst Case)
        Combined steps → O(n(n + m))
        The specific explanation is in the PDF.
        """
        # TODO
        n = G.vertex_num
        sd = [0] * n

        # Calculate the core number for each node.
        core_num = kCoreBaseStructuralDiversity.coreNumber(G)

        # Construct a new adjacency list for each node whose core number is greater than or equal to (k+1).
        adj_k1core = kCoreBaseStructuralDiversity.coreAdj(G, core_num, k + 1)

        for v in range(n):
            if core_num[v] < k + 1 or not adj_k1core[v]:
                sd[v] = 0
                continue

            neighbor_nodes = adj_k1core[v]
            #Construct the adjacency subgraph for each filtered node.
            sub_graph = kCoreBaseStructuralDiversity.buildSubgraph(neighbor_nodes, adj_k1core)
            #Calculate the number of connected components that meet the K-core requirements in the subgraph.
            sd[v] = kCoreBaseStructuralDiversity.kCoreCount(sub_graph, k)

        return sd


    # Construct a new adjacency list for each node whose core number is greater than or equal to (k+1), excluding the node itself.
    @staticmethod
    def coreAdj(G, core_num, threshold):
        n = G.vertex_num
        adj_k1core = [[] for _ in range(n)]
        adj_k1core = {u: [v for v in G.adj_list[u] if core_num[v] >= threshold]
              for u in range(n) if core_num[u] >= threshold}
        return adj_k1core
    #Construct the adjacency subgraph with this node removed.
    @staticmethod
    def buildSubgraph(neighbor_nodes, adj_k1core):
        neighbor_nodes = set(neighbor_nodes)
        sub_graph = {u: [] for u in neighbor_nodes}
        graph_edges = {(u, w) for u in neighbor_nodes for w in adj_k1core[u] if w in neighbor_nodes and u < w}
        for u, w in graph_edges:
            sub_graph[u].append(w)
            sub_graph[w].append(u)
        return sub_graph

    #Calculate the core number of each node using the core number decomposition algorithm based on a flat array.
    @staticmethod
    def initialize_data_structures(G, n):
        degree = list(map(len, G.adj_list))
        max_deg = max(degree, default=0)
        degree_fre = [0] * (max_deg + 1)
        initial_start = [0] * (max_deg + 1)
        vert_sorted = [0] * n
        vertex_current_positions = [0] * n
        core_num = degree[:]
        for i in range(len(degree)):
            degree_fre[degree[i]] += 1

        for i, prev_sum in enumerate(initial_start[:-1], 1):
            initial_start[i] = prev_sum + degree_fre[i - 1]

        current_begin = initial_start.copy()

        for u in range(n):
            vertex_current_positions[u], vert_sorted[current_begin[degree[u]]], current_begin[degree[u]] = \
                current_begin[degree[u]], u, current_begin[degree[u]] + 1

        return degree, max_deg, degree_fre, initial_start, vert_sorted, vertex_current_positions, current_begin, core_num
    @staticmethod
    def coreNumber(G):
        n = G.vertex_num
        if n == 0: return []

        degree, max_deg, degree_fre, initial_start, vert_sorted, vertex_current_positions, current_begin, core_num = \
        kCoreBaseStructuralDiversity.initialize_data_structures(G, n)
        core_num = degree[:]
        for i in range(n):
            u, cur_core = vert_sorted[i], core_num[vert_sorted[i]]
            for v in G.adj_list[u]:
                if core_num[v] >  cur_core :
                    target_pos = initial_start[core_num[v]]
                    swap_point = vert_sorted[target_pos]
                    neigh_pos = vertex_current_positions[v]
                    if v != swap_point: vert_sorted[neigh_pos], vert_sorted[target_pos], vertex_current_positions[v], vertex_current_positions[swap_point] = vert_sorted[target_pos], vert_sorted[neigh_pos], target_pos, neigh_pos
                    initial_start[core_num[v]] += 1
                    core_num[v] -= 1

        return core_num



    #Select the connected components that satisfy the k-core condition in the subgraph of the undirected graph.
    @staticmethod
    def _init_kcore_variable(sub_graph, k):
        node_deg = {node: len(neigh) for node, neigh in sub_graph.items()}
        node_pre = [node for node, deg in node_deg.items() if deg < k]
        del_node = set()
        return node_deg, node_pre, del_node
    @staticmethod
    def kCoreCount(sub_graph, k):
        if not sub_graph: return 0

        node_deg, node_pre, del_node = kCoreBaseStructuralDiversity._init_kcore_variable(sub_graph, k)

        idx = 0
        while idx < len(node_pre):
            cur_node = node_pre[idx]
            idx += 1
            if cur_node not in del_node: del_node.add(cur_node)
            else: continue
            del_node.add(cur_node)
            for neigh in sub_graph[cur_node]:
                if neigh in del_node: continue
                node_deg[neigh] -= 1
                if node_deg[neigh] == k - 1:
                    node_pre.append(neigh)

        remain_node = list(filter(lambda cur_node: cur_node not in del_node, sub_graph))
        if not remain_node: return 0

        return kCoreBaseStructuralDiversity.ccCount(sub_graph, remain_node)

    #Calculate the number of connected components.
    @staticmethod
    def ccCount(sub_graph, remain_node):
        remain_set = set(remain_node)
        visited = set()
        component_count = 0

        def dfs(node):
            visited.add(node)
            neighbors = sub_graph.get(node, [])
            for neigh in neighbors:
                if neigh in remain_set and neigh not in visited:
                    dfs(neigh)

        for cur_node in remain_set:
            if cur_node not in visited:
                component_count += 1
                dfs(cur_node)
        return component_count


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
