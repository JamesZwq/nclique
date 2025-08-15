#!/usr/bin/env python3
# Auto-generated for 5582348

STUDENT_ID = "5582348"
STUDENT_NAME = "Wanyun Wu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n  # τ_k(v) for each vertex

        for v in range(n):
            neighbors = set(G.adj_list[v])
            if not neighbors:
                continue

            sub_edges = []
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbors and u < w:
                        sub_edges.append((u, w))

            sub_adj = {u: set() for u in neighbors}
            for u, w in sub_edges:
                sub_adj[u].add(w)
                sub_adj[w].add(u)

            k_core_nodes = kCoreBaseStructuralDiversity._extract_k_core(sub_adj, k)

            sd[v] = kCoreBaseStructuralDiversity._count_components(sub_adj, k_core_nodes)

        return sd

    @staticmethod
    def _extract_k_core(adj, k):
        deg = {u: len(adj[u]) for u in adj}
        q = deque([u for u in adj if deg[u] < k])
        removed = set()

        while q:
            u = q.popleft()
            if u in removed:
                continue
            removed.add(u)
            for v in adj[u]:
                if v not in removed:
                    deg[v] -= 1
                    if deg[v] == k - 1:
                        q.append(v)

        return {u for u in adj if u not in removed}

    @staticmethod
    def _count_components(adj, nodes):
        visited = set()
        count = 0

        for u in nodes:
            if u not in visited:
                count += 1
                q = deque([u])
                visited.add(u)
                while q:
                    x = q.popleft()
                    for y in adj[x]:
                        if y in nodes and y not in visited:
                            visited.add(y)
                            q.append(y)
        return count

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
