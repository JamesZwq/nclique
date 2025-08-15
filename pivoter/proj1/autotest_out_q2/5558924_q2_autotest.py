#!/usr/bin/env python3
# Auto-generated for 5558924

STUDENT_ID = "5558924"
STUDENT_NAME = "Selina Zhang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from collections import defaultdict
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

        for v in range(n):
            # Process each vertex
            sd[v] = kCoreBaseStructuralDiversity._process_single_vertex(G, v, k)
        return sd
    @staticmethod
    def _process_single_vertex(G, v, k):
        adj = G.adj_list[v]# Get neighbors of vertex v
        if k <= 0:
            return 0   # If k is invalid, return 0
        if len(adj) == 0:
            return 0   # If v has no neighbors, diversity is 0

        induced_subgraph = defaultdict(list) # Adjacency list for induced subgraph of neighbors
        adj_set = set(adj)

        for u in adj_set:
            # Add edges among neighbors to the induced subgraph
            induced_subgraph[u].extend(
                [w for w in G.adj_list[u] if w in adj_set]
            )
        return kCoreBaseStructuralDiversity._count_k_core_components(induced_subgraph, k)

    @staticmethod
    def _count_k_core_components(induced_subgraph, k):
        if not induced_subgraph:
            return 0 # If subgraph is empty, return 0
        # Filter out nodes with degree less than k
        active = kCoreBaseStructuralDiversity._filter_k_core_nodes(induced_subgraph, k)
        if not active:
            return 0
        # Count the number of connected components in the remaining graph
        return kCoreBaseStructuralDiversity._count_connected_components(induced_subgraph, active)


    @staticmethod
    def _filter_k_core_nodes(induced_subgraph, k):
        node_degrees = {u: len(induced_subgraph[u]) for u in induced_subgraph} # Initial degree of each node
        queue = deque([u for u in induced_subgraph if node_degrees[u] < k]) # Nodes to be removed
        active_nodes = set(induced_subgraph.keys()) - set(queue)  # Initial k-core candidates

        while queue:
            u = queue.popleft()
            for v in induced_subgraph[u]:
                if v in active_nodes:
                    node_degrees[v] -= 1
                    if node_degrees[v] < k:
                        active_nodes.remove(v)
                        queue.append(v)

        return active_nodes # Final set of nodes in the k-core

    @staticmethod
    def _count_connected_components(induced_subgraph, node_set):
        discovered = set()
        total = 0
        pending = list(node_set)

        while pending:
            start = pending.pop()
            if start in discovered:
                continue

            total += 1 # New connected component found
            stack = [start]

            while stack:
                node = stack.pop()
                if node in discovered:
                    continue
                discovered.add(node)

                neighbors =induced_subgraph.get(node, [])
                for nxt in neighbors:
                    if nxt in node_set and nxt not in discovered:
                        stack.append(nxt)

        return total # Number of connected components


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
