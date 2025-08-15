#!/usr/bin/env python3
# Auto-generated for 5495326

STUDENT_ID = "5495326"
STUDENT_NAME = "Yichen Wang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity:
    @staticmethod
    def process(G, k):
        n, adj = G.vertex_num, G.adj_list
        sd = [0]*n

        for v in range(n):
            neigh = adj[v]
            d = len(neigh)
            if d < k:
                continue
            Nv = set(neigh)

            # build induced subgraph
            sub = {u: [] for u in neigh}
            for u in neigh:
                for w in adj[u]:
                    if w in Nv:
                        sub[u].append(w)

            # compute k-core
            degrees = {u: len(sub[u]) for u in sub}
            deleted = set()
            todo = deque(u for u,d0 in degrees.items() if d0 < k)
            while todo:
                u = todo.popleft()
                if u in deleted:
                    continue
                deleted.add(u)
                for w in sub[u]:
                    if w not in deleted:
                        degrees[w] -= 1
                        if degrees[w] == k-1:
                            todo.append(w)

            # collect core nodes
            core = [u for u in sub if u not in deleted]
            if not core:
                continue

            # count components
            visited, counts = set(), 0
            for u in core:
                if u not in visited:
                    counts += 1
                    dq = deque([u])
                    visited.add(u)
                    while dq:
                        x = dq.popleft()
                        for w in sub[x]:
                            if w not in deleted and w not in visited:
                                visited.add(w)
                                dq.append(w)
            sd[v] = counts

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
