#!/usr/bin/env python3
# Auto-generated for 5522276

STUDENT_ID = "5522276"
STUDENT_NAME = "Junyi Zheng"

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
        # Iterate over all v_node and calculate the number of k-cores
        # of each node's neighboring subgraphs separately
        # Check if this vertex has any neighbors
        # sd[v] is the number of subgraphs that satisfy the k-core condition,
        # e.g., the degree of every point inside this subgraph should be greater than or equal to k
        for v in range(n):
            # Get the neighbor of vertex v and form a subgraph
            node_neighbours = G.adj_list[v]
            if not node_neighbours:
                sd[v] = 0
                continue

            # subgraph does not contain v_node v itself is handled with adjacency list
            child_G = kCoreBaseStructuralDiversity.build_child_G(G, node_neighbours)
            # Call k_core_count_number to compute the number of k-cores in a sub graph
            sd[v] = kCoreBaseStructuralDiversity.k_core_count_number(child_G, k)
        return sd

    @staticmethod
    def k_core_count_number(child_G, k):
        # Delete any point in the subgraph
        # if it does not satisfy the condition greater than or equal to k.
        # Delete points with unsatisfactory node_degrees
        if not child_G:
            return 0
        after_split_deal_G, delete = kCoreBaseStructuralDiversity.k_core_nodes_calculate(child_G, k)
        after_split_deal_G = []
        for u in child_G:
            if u not in delete:
                after_split_deal_G.append(u)
        if not after_split_deal_G:
            return 0

        visited = set()
        k_core_count = 0
        # DFS
        for node in after_split_deal_G:
            if node not in visited:
                k_core_count += 1
                stack = [node]
                visited.add(node)
                # dfs
                while stack:
                    u = stack.pop()
                    for v in child_G[u]:
                        if v not in delete and v not in visited:
                            visited.add(v)
                            stack.append(v)

        return k_core_count

    @staticmethod
    def build_child_G(G, node_neighbours):
        # Construct induced subgraph (only edges
        # between node_neighbours, not original v_node v)
        # de-duplicate
        node_neighbours = set(node_neighbours)
        child_G = {}
        for u in node_neighbours:
            child_G[u] = []
        for u in node_neighbours:
            for u_neighbours in G.adj_list[u]:
                # Avoid joining non-neighboring points
                if u_neighbours in node_neighbours:
                    child_G[u].append(u_neighbours)
        return child_G
    @staticmethod
    def k_core_nodes_calculate(child_G, k):
        # node_degrees = {u: len(node_neighbours) for u, node_neighbours in child_G.items()}
        node_degrees = {}
        for u in child_G:
            node_degrees[u] = len(child_G[u])
        delete = set()
        # Find all v_node with degree less than k and process them,
        # since they do not satisfy the requirement of being greater than or equal to k-core
        v_node = deque()
        for u in child_G:
            if node_degrees[u] < k:
                v_node.append(u)
                delete.add(u)
        # After the unsatisfied v_node are processed,
        # the degree of the remaining v_node will also change and the degree needs to be updated.
        while v_node:
            u = v_node.popleft()
            for v in child_G[u]:
                if v not in delete:
                    node_degrees[v] -= 1
                    # Determine if the degree of the remaining nodes <= k If less than it is also processed.
                    if node_degrees[v] < k:
                        delete.add(v)
                        v_node.append(v)

        after_split_deal_G = []
        for u in child_G:
            if u not in delete:
                after_split_deal_G.append(u)
        return after_split_deal_G, delete


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
