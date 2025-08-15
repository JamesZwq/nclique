#!/usr/bin/env python3
# Auto-generated for 5504446

STUDENT_ID = "5504446"
STUDENT_NAME = "Zhibo Zhang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _build_induced_subgraph(G, nodes):
        """
        Constructs the adjacency list of the subgraph induced by nodes.
        Returns: {node: set(neighbors)}
        """
        sub_adj = {u: set() for u in nodes}  # Initialize the subgraph adjacency list
        for u in nodes:
            for v in G.adj_list[u]:
                if v in nodes:
                    sub_adj[u].add(v)  # Only keep edges between neighbors
        return sub_adj

    @staticmethod
    def _k_core_peeling(adj, k):
        """
        Recursively remove nodes with degree less than k from the adjacency list adj.
        Return the removed adjacency list
        """
        changed = True
        while changed:
            changed = False
            # find all nodes with degree less than k
            to_remove = [u for u in adj if len(adj[u]) < k]
            if to_remove:
                changed = True
                for u in to_remove:
                    for v in adj[u]:
                        adj[v].discard(u)
                    del adj[u]
        # adj only contains nodes with degree >= k after stripping
        return adj

    @staticmethod
    def _count_connected_components(adj):
        """
        Count the number of connected components in the adjacency list adj
        """
        visited = set()
        count = 0
        for u in adj:
            if u not in visited:
                # BFS traversal finds a connected component
                queue = deque([u])
                visited.add(u)
                while queue:
                    curr = queue.popleft()
                    for v in adj[curr]:
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)
                count += 1  # Each time a connected component is found, the count is increased by one
        return count  

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
        sd = [0] * n  # Stores the k-core structure diversity of each vertex
        for v in range(n):
            neighbors = G.adj_list[v]  # Get the neighbor set of vertex v
            if not neighbors:
                sd[v] = 0  # No neighbors, structural diversity is 0
                continue
            # Constructing neighbor-induced subgraphs
            sub_nodes = set(neighbors) # Neighbor node set
            sub_adj = kCoreBaseStructuralDiversity._build_induced_subgraph(G, sub_nodes)
            # k-core stripping
            adj = {u: set(sub_adj[u]) for u in sub_adj}  # make a copy for stripping
            kCoreBaseStructuralDiversity._k_core_peeling(adj, k)
            # Count the remaining connected components
            count = kCoreBaseStructuralDiversity._count_connected_components(adj)
            sd[v] = count  # Record the k-core structure diversity of vertex v
        return sd

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
