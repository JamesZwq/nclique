#!/usr/bin/env python3
# Auto-generated for 5505424

STUDENT_ID = "5505424"
STUDENT_NAME = "Xinran Yang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # Function to find connected components in a given subgraph
    @staticmethod
    def find_connected_components(subgraph_nodes, subgraph_adj):
        visited_subgraph = {node: False for node in subgraph_nodes}
        component_count = 0

        for start_node in subgraph_nodes:
            if not visited_subgraph[start_node]:
                component_count += 1
                q = deque([start_node])
                visited_subgraph[start_node] = True
                while q:
                    u = q.popleft()
                    for v_neighbor in subgraph_adj.get(u, []):
                        if v_neighbor in subgraph_nodes and not visited_subgraph[v_neighbor]:
                            visited_subgraph[v_neighbor] = True
                            q.append(v_neighbor)
        return component_count

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
        tau = [0] * n
        if n == 0:
            return tau

        adj = G.adj_list

        # Main loop over all vertices
        for v_idx in range(n):
            N_v = set(adj[v_idx])  # Get neighbors of current vertex v_idx

            if not N_v:
                tau[v_idx] = 0
                continue

            # Construct the neighbor-induced subgraph
            neighbor_induced_adj = {node: [] for node in N_v}
            for u in N_v:
                for neighbor_of_u in adj[u]:
                    if neighbor_of_u in N_v:
                        neighbor_induced_adj[u].append(neighbor_of_u)

            # For k=0, it's just the number of connected components in the neighbor-induced subgraph
            if k == 0:
                tau[v_idx] = kCoreBaseStructuralDiversity.find_connected_components(list(N_v), neighbor_induced_adj)
                continue

            # For k > 0, perform k-core decomposition
            current_nodes_in_core = set(N_v)

            # Calculate initial degrees within the neighbor-induced subgraph
            current_degrees = {node: len(neighbor_induced_adj[node]) for node in current_nodes_in_core}

            q_peeling = deque()
            for node in current_nodes_in_core:
                if current_degrees[node] < k:
                    q_peeling.append(node)

            while len(q_peeling) > 0:
                node_to_remove = q_peeling.popleft()
                # Ensure it hasn't been removed by another path
                if node_to_remove not in current_nodes_in_core:
                    continue
                current_nodes_in_core.remove(node_to_remove)
                for neighbor_of_removed in neighbor_induced_adj.get(node_to_remove, []):
                    if neighbor_of_removed in current_nodes_in_core:
                        current_degrees[neighbor_of_removed] -= 1
                        if current_degrees[neighbor_of_removed] < k:
                            q_peeling.append(neighbor_of_removed)

            # Reconstruct adjacency list for the k-core (remaining nodes)
            k_core_adj = {node: [] for node in current_nodes_in_core}
            for u in current_nodes_in_core:
                for neighbor_of_u in neighbor_induced_adj.get(u, []):
                    if neighbor_of_u in current_nodes_in_core:
                        k_core_adj[u].append(neighbor_of_u)

            # Count connected components in the k-core
            tau[v_idx] = kCoreBaseStructuralDiversity.find_connected_components(list(current_nodes_in_core), k_core_adj)

        return tau


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
