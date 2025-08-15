#!/usr/bin/env python3
# Auto-generated for 5528317

STUDENT_ID = "5528317"
STUDENT_NAME = "Yankun Wei"

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
        sd = [0] * n

        for node in range(n):
            # Compute the neighbour-induced subgraph
            sub_graph = kCoreBaseStructuralDiversity.__cac_nbr_graph(G, node)
            if sub_graph is None:
                continue

            # Calculate k-core-based structural diversity
            k_cores_nodes = kCoreBaseStructuralDiversity.__calc_k_cores(sub_graph, k)
            sd_value = kCoreBaseStructuralDiversity.__find_nb_of_connected_components(sub_graph, k_cores_nodes)
            sd[node] = sd_value
        
        return sd


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def __cac_nbr_graph(G, node):
        """Calculate the neighbour-induced subgraph"""
        # Find its neighbors
        neighbors = set(G.adj_list[node])
        if not neighbors:
            return None

        # Compute the neighbour-induced subgraph
        sub_graph = {}
        for u in neighbors:
            sub_graph[u] = []
            for v in G.adj_list[u]:
                if v not in neighbors:
                    continue
                # Add one edge into the subgraph
                sub_graph[u].append(v)
        
        return sub_graph
        

    @staticmethod
    def __calc_k_cores(g, k):
        """Calculate the k-core of the graph"""
        # Compute the degrees of all nodes
        degrees = {node: len(g[node]) for node in g}

        # Find all nodes with degree less than k.
        excluded_nodes = set()
        q = deque()
        for node in g:
            if degrees[node] >= k:
                continue
            excluded_nodes.add(node)
            q.append(node)

        while q:
            node = q.popleft()
            for neighbor in g[node]:
                if neighbor in excluded_nodes:
                    continue
                # The neighbor's degree decreases by 1
                degrees[neighbor] = degrees[neighbor] - 1
                if degrees[neighbor] >= k:
                    continue
                # The node is not in any k-core
                excluded_nodes.add(neighbor)
                q.append(neighbor)

        # The remaining nodes are in k-core
        graph_nodes = set(g.keys())
        k_cores_nodes = graph_nodes - excluded_nodes
        return k_cores_nodes
    
    @staticmethod
    def __find_nb_of_connected_components(g, nodes):
        """Find the number of connected components"""
        visited_nodes = set()
        def bfs(start_node):
            """Breadth-first search"""
            visited_nodes.add(start_node)
            q = deque([start_node])
            cc_nodes = set([start_node])
            while q:
                node = q.popleft()
                # Visit its neighbor
                for neighbor in g[node]:
                    if neighbor not in nodes:
                        continue
                    if neighbor in visited_nodes:
                        continue
                    visited_nodes.add(neighbor)
                    cc_nodes.add(neighbor)
                    q.append(neighbor)

            return cc_nodes

        nb_cc = 0
        for node in nodes:
            if node in visited_nodes:
                continue

            # Use breadth-first search to find one connected component
            cc_nodes = bfs(node)
            nb_cc = nb_cc + 1

        return nb_cc

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
