#!/usr/bin/env python3
# Auto-generated for 5459829

STUDENT_ID = "5459829"
STUDENT_NAME = "Cheng Li"

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
        n_nodes = G.vertex_num
        diversity_counts = [0] * n_nodes
        for idx in range(n_nodes):
            neighbors = set(G.adj_list[idx])
            if not neighbors:
                diversity_counts[idx] = 0
                continue
            # Induced subgraph adjacency construction (bidirectional for undirected)
            induced_adj = {node: [] for node in neighbors}
            for node in neighbors:
                for nbr in G.adj_list[node]:
                    if nbr in neighbors and node != nbr:
                        if nbr not in induced_adj[node]:
                            induced_adj[node].append(nbr)
                        if node not in induced_adj[nbr]:
                            induced_adj[nbr].append(node)
            diversity_counts[idx] = kCoreBaseStructuralDiversity._num_kcore_components(induced_adj, k)
        return diversity_counts

    @staticmethod
    def _num_kcore_components(graph, k):
        # Compute all k-core vertices (prune low-degree iteratively)
        if not graph:
            return 0
        deg = {u: len(adj) for u, adj in graph.items()}
        removed = set()
        to_prune = deque(u for u in graph if deg[u] < k)
        removed.update(to_prune)
        while to_prune:
            node = to_prune.popleft()
            for nbr in graph[node]:
                if nbr not in removed:
                    deg[nbr] -= 1
                    if deg[nbr] < k:
                        to_prune.append(nbr)
                        removed.add(nbr)
        # Remaining nodes form k-core(s)
        core_nodes = [u for u in graph if u not in removed]
        if not core_nodes:
            return 0
        # Count number of connected components in the k-core
        explored = set()
        component_total = 0
        for v in core_nodes:
            if v in explored:
                continue
            component_total += 1
            stack = deque([v])
            explored.add(v)
            while stack:
                curr = stack.popleft()
                for adj in graph[curr]:
                    if adj not in removed and adj not in explored:
                        explored.add(adj)
                        stack.append(adj)
        return component_total

    ################################################################################
    # Any extra helper functions should be defined inside this class.
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
