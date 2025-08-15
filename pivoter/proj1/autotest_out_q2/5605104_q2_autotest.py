#!/usr/bin/env python3
# Auto-generated for 5605104

STUDENT_ID = "5605104"
STUDENT_NAME = "Zhaoyu Wang"

# ======= 学生代码 =======
from collections import deque, defaultdict

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
        sd = [0] * n  # τ_k(v) for each vertex v

        # Build adjacency list from G
        adj = G.adj_list




        # adj = [[] for _ in range(n)]
        # for u, v in G.edges:
        #     adj[u].append(v)
        #     adj[v].append(u)

        def compute_k_cores(sub_adj):
            """Return a list of sets of nodes representing each k-core."""
            degree = {u: len(neigh) for u, neigh in sub_adj.items()}
            core = dict()
            bin = deque([u for u in degree if degree[u] < k])

            while bin:
                u = bin.popleft()
                for v in sub_adj[u]:
                    if v in degree:
                        degree[v] -= 1
                        if degree[v] < k and v not in bin:
                            bin.append(v)
                del degree[u]

            # Extract connected components from remaining nodes
            visited = set()
            components = []

            for u in degree:
                if u not in visited:
                    comp = set()
                    dq = deque([u])
                    visited.add(u)
                    while dq:
                        curr = dq.popleft()
                        comp.add(curr)
                        for nei in sub_adj[curr]:
                            if nei in degree and nei not in visited:
                                visited.add(nei)
                                dq.append(nei)
                    components.append(comp)

            return components

        # Process each vertex
        for v in range(n):
            neighbors = set()
            for u in adj[v]:
                neighbors.add(u)

            # Build neighbor-induced subgraph
            sub_adj = defaultdict(list)
            for u in neighbors:
                for w in adj[u]:
                    if w in neighbors:
                        sub_adj[u].append(w)

            # Get k-cores in neighbor-induced subgraph
            k_cores = compute_k_cores(sub_adj)
            sd[v] = len(k_cores)

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
