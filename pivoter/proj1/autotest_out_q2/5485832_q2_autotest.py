#!/usr/bin/env python3
# Auto-generated for 5485832

STUDENT_ID = "5485832"
STUDENT_NAME = "Shengze Kou"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque

################################################################################
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    def process(G, k):
        # Count number of k-core components in a given adjacency list
        def get_k_core_count(adj):
            visited = set()
            components = []

            # Extract connected components using DFS
            def dfs(u, comp):
                visited.add(u)
                comp.append(u)
                for v in adj[u]:
                    if v not in visited:
                        dfs(v, comp)

            for node in adj:
                if node not in visited:
                    comp = []
                    dfs(node, comp)
                    components.append(comp)

            count = 0
            for comp in components:
                # Build subgraph for each connected component
                sub_adj = {u: [v for v in adj[u] if v in comp] for u in comp}

                # Iteratively remove nodes with degree < k
                changed = True
                while changed:
                    changed = False
                    remove = [u for u in sub_adj if len(sub_adj[u]) < k]
                    if remove:
                        changed = True
                        for u in remove:
                            for v in sub_adj[u]:
                                sub_adj[v].remove(u)
                            del sub_adj[u]

                # Count remaining k-core connected components
                if sub_adj:
                    visited2 = set()
                    for node in sub_adj:
                        if node not in visited2:
                            stack = [node]
                            visited2.add(node)
                            while stack:
                                u = stack.pop()
                                for v in sub_adj[u]:
                                    if v not in visited2:
                                        visited2.add(v)
                                        stack.append(v)
                            count += 1
            return count

        n = G.vertex_num
        result = [0] * n  # τ_k(v) for each vertex

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                continue

            # Build neighbor-induced subgraph
            sub_nodes = set(neighbors)
            sub_adj = {u: [] for u in sub_nodes}
            for u in sub_nodes:
                for w in G.adj_list[u]:
                    if w in sub_nodes:
                        sub_adj[u].append(w)

            # Count k-core components in subgraph
            result[v] = get_k_core_count(sub_adj)

        return result

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
