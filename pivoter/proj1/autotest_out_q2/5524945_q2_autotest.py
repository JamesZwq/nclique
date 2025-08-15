#!/usr/bin/env python3
# Auto-generated for 5524945

STUDENT_ID = "5524945"
STUDENT_NAME = "Adhishwar Narayan Tiwari"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _find_k_core_components(num_nodes, adj_list, k):
        """
        Finding how many k-core groups are in a graph.

        Starting by removing nodes that do not have enough connections.
        Then counting how many separate groups are left.

        Parameters:
        - num_nodes (int): How many nodes we have
        - adj_list (List[List[int]]): List showing which nodes connect to which
        - k (int): The minimum connections needed

        Returns:
        - int: Number of k-core groups found
        """
        # Handling special case when k=0
        # When k=0, just counting all connected groups
        if k == 0:
            if not num_nodes:
                return 0
            visited = [False] * num_nodes
            components = 0

            for i in range(num_nodes):
                if not visited[i]:
                    components += 1
                    q = deque([i])
                    visited[i] = True

                    while q:
                        u = q.popleft()
                        for v_neighbor in adj_list[u]:
                            if not visited[v_neighbor]:
                                visited[v_neighbor] = True
                                q.append(v_neighbor)
            return components

        # Checking if k-core is possible
        # Need more nodes than k to have any k-cores
        if num_nodes <= k:
            return 0

        # Starting k-core decomposition, removing weak nodes
        degrees = [len(adj) for adj in adj_list]
        removed = [False] * num_nodes
        q = deque()

        # Finding nodes with too few connections
        for i in range(num_nodes):
            if degrees[i] < k:
                q.append(i)
                removed[i] = True

        # Removing nodes and updating their neighbors
        while q:
            u = q.popleft()
            for v_neighbor in adj_list[u]:
                if not removed[v_neighbor]:
                    degrees[v_neighbor] -= 1
                    if degrees[v_neighbor] < k:
                        q.append(v_neighbor)
                        removed[v_neighbor] = True

        # Counting connected groups in remaining nodes
        components = 0
        visited = [False] * num_nodes

        for i in range(num_nodes):
            # Checking if node is still there and not visited
            if not removed[i] and not visited[i]:
                components += 1
                # Starting search to find all nodes in this group
                component_q = deque([i])
                visited[i] = True

                while component_q:
                    u_comp = component_q.popleft()
                    # Looking at neighbors that are also still there
                    for v_comp_neighbor in adj_list[u_comp]:
                        if not removed[v_comp_neighbor] and not visited[v_comp_neighbor]:
                            visited[v_comp_neighbor] = True
                            component_q.append(v_comp_neighbor)

        return components

    @staticmethod
    def process(G, k):
        """
        Computing structural diversity for every node in the graph.

        Parameters:
        G : UndirectedUnweightedGraph
            The graph we're working with
        k : int
            The minimum connections needed for k-core

        Returns:
        List[int] = List where each position shows the diversity value for that node
        """
        n = G.vertex_num
        tau_values = [0] * n

        # Going through each node to calculate its diversity
        for v_idx in range(n):
            neighbors = G.adj_list[v_idx]
            num_neighbors = len(neighbors)

            # Checking if node has enough neighbors
            # If not enough neighbors, no k-core can exist
            if k > 0 and num_neighbors <= k:
                tau_values[v_idx] = 0
                continue

            # Handling nodes with no neighbors
            if not neighbors:
                tau_values[v_idx] = 0
                continue

            # Building subgraph from neighbors only
            # Making new IDs for neighbors like 0,1,2....
            node_map = {node_id: i for i, node_id in enumerate(neighbors)}
            subgraph_adj_list = [[] for _ in range(num_neighbors)]

            # Using set for fast checking if node is a neighbor
            neighbor_set = set(neighbors)

            for i, u in enumerate(neighbors):
                # Looking at each neighbor's connections
                for w in G.adj_list[u]:
                    # Checking if this connection is also our neighbor
                    if w in neighbor_set:
                        j = node_map[w]
                        subgraph_adj_list[i].append(j)

            # Finding k-core groups in the neighbor subgraph
            num_components = kCoreBaseStructuralDiversity._find_k_core_components(
                num_neighbors, subgraph_adj_list, k
            )
            tau_values[v_idx] = num_components

        return tau_values


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
