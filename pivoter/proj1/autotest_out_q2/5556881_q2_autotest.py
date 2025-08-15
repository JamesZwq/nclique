#!/usr/bin/env python3
# Auto-generated for 5556881

STUDENT_ID = "5556881"
STUDENT_NAME = "(Leo) Lixian Deng"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):

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
        adj = G.adj_list
        τ = [0] * n

        for v in range(n):
            neighbors = set(adj[v])
            if not neighbors:
                continue

            # Step 1: Build neighbor-induced subgraph (excluding v)
            sub_adj = kCoreBaseStructuralDiversity._build_subgraph(adj, neighbors)

            # Step 2: Prune to k-core using degree peeling
            kcore_nodes = kCoreBaseStructuralDiversity._extract_k_core(sub_adj, k)

            # Step 3: Count connected components using DFS
            τ[v] = kCoreBaseStructuralDiversity._count_connected_components(sub_adj, kcore_nodes)

        return τ

    @staticmethod
    def _build_subgraph(adj, node_set):
        sub_adj = {u: [] for u in node_set}
        for u in node_set:
            for v in adj[u]:
                if v in node_set:
                    sub_adj[u].append(v)
        return sub_adj

    @staticmethod
    def _extract_k_core(sub_adj, k):
        deg = {u: len(sub_adj[u]) for u in sub_adj}
        queue = deque([u for u in deg if deg[u] < k])
        removed = set()

        while queue:
            u = queue.popleft()
            removed.add(u)
            for nbr in sub_adj[u]:
                if nbr not in removed and nbr in deg:
                    deg[nbr] -= 1
                    if deg[nbr] == k - 1:
                        queue.append(nbr)
            del deg[u]

        return set(deg.keys())

    @staticmethod
    def _count_connected_components(sub_adj, valid_nodes):
        visited = set()
        count = 0
        for u in valid_nodes:
            if u in visited:
                continue
            count += 1
            stack = [u]
            visited.add(u)
            while stack:
                curr = stack.pop()
                for nbr in sub_adj[curr]:
                    if nbr in valid_nodes and nbr not in visited:
                        visited.add(nbr)
                        stack.append(nbr)
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
