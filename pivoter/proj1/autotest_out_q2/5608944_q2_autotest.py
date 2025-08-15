#!/usr/bin/env python3
# Auto-generated for 5608944

STUDENT_ID = "5608944"
STUDENT_NAME = "Tong Zhang"

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
        tau_values = [0 for _ in range(n)]  # Store k-core diversity value for each vertex

        for idx in range(n):
            neighbors = G.adj_list[idx]

            # Skip isolated vertex (no neighbors)
            if len(neighbors) == 0:
                tau_values[idx] = 0
                continue

            # Construct a subgraph induced by the neighbors of vertex idx
            local = set(neighbors)  # Set of neighbors (excluding idx itself)
            subgraph = dict()

            # Initialize adjacency list for the subgraph
            for u in local:
                subgraph[u] = []

            # Populate subgraph with edges between neighbors
            for u in local:
                for v in G.adj_list[u]:
                    if v in local:
                        subgraph[u].append(v)

            # Compute number of connected components in k-core of the subgraph
            tau_values[idx] = kCoreBaseStructuralDiversity._kcore_connected_count(subgraph, k)

        return tau_values

    @staticmethod
    def _kcore_connected_count(graph, k):
        """
        Compute the number of connected components in the k-core of a given graph.
        A node is in the k-core if its degree is >= k after recursive pruning.
        """
        if len(graph) == 0:
            return 0

        #  Prune all nodes whose degree is initially less than k
        deg = dict((u, len(v)) for u, v in graph.items())  # Degree map
        removed_nodes = set()
        pruning = deque(node for node in graph if deg[node] < k)
        removed_nodes.update(pruning)

        # Recursively remove nodes and update neighbors' degrees
        while pruning:
            u = pruning.popleft()
            for v in graph[u]:
                if v not in removed_nodes:
                    deg[v] -= 1  # Decrease degree
                    if deg[v] < k:
                        removed_nodes.add(v)
                        pruning.append(v)

        # Remaining nodes after pruning are considered part of the k-core
        residual_nodes = [u for u in graph if u not in removed_nodes]
        if len(residual_nodes) == 0:
            return 0

        # Count connected components in the remaining graph using DFS
        visited = set()
        component_total = 0

        for root in residual_nodes:
            if root in visited:
                continue

            # New connected component found
            component_total += 1
            stack = list([root])
            visited.add(root)

            # Perform DFS to visit all nodes in the current component
            while stack:
                current = stack.pop()
                for neighbor in graph[current]:
                    if neighbor not in removed_nodes and neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)

        return component_total

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
