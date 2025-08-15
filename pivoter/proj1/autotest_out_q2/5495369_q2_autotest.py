#!/usr/bin/env python3
# Auto-generated for 5495369

STUDENT_ID = "5495369"
STUDENT_NAME = "Chloe Wang"

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
        # get neighbors for v
        for v in range(n):
            nbr = set(G.adj_list[v])
            if not nbr:
                sd[v] = 0
                continue

            neighbors = list(nbr)
            node_index = {u: i for i, u in enumerate(neighbors)}
            subgraph_edges = []
            # loop all edges
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in nbr and node_index[u] < node_index[w]:
                        subgraph_edges.append((node_index[u], node_index[w]))
            # return 0 if no subgraph
            if not subgraph_edges:
                sd[v] = 0
                continue
            # generate Neighbour-induced subgraph
            H = UndirectedUnweightedGraph([(len(neighbors), len(subgraph_edges))] + subgraph_edges)
            # get k_core for subgraph
            k_core_nodes = kCoreBaseStructuralDiversity.generate_k_core(H, k)
            # count conected_component
            sd[v] = kCoreBaseStructuralDiversity.get_CC(H, k_core_nodes)

        return sd

    @staticmethod
    def generate_k_core(G, k):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        Set[int]: Set of node indices that belong to the k-core of the graph
        """
        n = G.vertex_num
        # get degree for all v
        degree = [len(G.adj_list[v]) for v in range(n)]
        valid_nodes = [True] * n
         # Iteratively remove nodes with degree < k
        changed = True
        while changed:
            changed = False
            for v in range(n):
                if valid_nodes[v] and degree[v] < k:
                    valid_nodes[v] = False
                    changed = True
                    for u in G.adj_list[v]:
                        if valid_nodes[u]:
                            degree[u] -= 1
        # return nodes are good
        return {v for v in range(n) if valid_nodes[v]}

    @staticmethod
    def get_CC(G, valid_nodes):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        valid_nodes: Set[int]
        Returns
        -------
        int: Number of connected components among valid nodes
        """
        visited = set()
        count = 0
        for v in valid_nodes:
            if v not in visited:
                queue = [v]
                visited.add(v)
                # BFS to traverse current connected_component
                while queue:
                    u = queue.pop()
                    for nei in G.adj_list[u]:
                        if nei in valid_nodes and nei not in visited:
                            visited.add(nei)
                            queue.append(nei)
                count += 1
        return count



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
