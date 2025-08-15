#!/usr/bin/env python3
# Auto-generated for 5519319

STUDENT_ID = "5519319"
STUDENT_NAME = "Yi Li"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n
        if n == 0:
            return sd
        inS_global = [False] * n
        deg_H_global = [0] * n

        for v in range(n):
            neighbors_v = G.adj_list[v]
            deg_v = len(neighbors_v)
            if deg_v == 0:
                sd[v] = 0
                continue
            for u in neighbors_v:
                inS_global[u] = True
            for u in neighbors_v:
                count = 0
                for w in G.adj_list[u]:
                    if inS_global[w]:
                        count += 1
                deg_H_global[u] = count
            q = deque()
            for u in neighbors_v:
                if deg_H_global[u] < k:
                    q.append(u)
            while q:
                u = q.popleft()
                if not inS_global[u]:
                    continue
                inS_global[u] = False

                for w in G.adj_list[u]:
                    if inS_global[w]:
                        deg_H_global[w] -= 1
                        if deg_H_global[w] < k:
                            q.append(w)
            core_set = set()
            for u in neighbors_v:
                if inS_global[u]:
                    core_set.add(u)
            visited_set = set()
            comp_count = 0

            for u in core_set:
                if u not in visited_set:
                    comp_count += 1
                    queue_comp = deque([u])
                    visited_set.add(u)

                    while queue_comp:
                        node = queue_comp.popleft()
                        for w in G.adj_list[node]:
                            if w in core_set and w not in visited_set:
                                visited_set.add(w)
                                queue_comp.append(w)

            sd[v] = comp_count
            for u in neighbors_v:
                inS_global[u] = False

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
