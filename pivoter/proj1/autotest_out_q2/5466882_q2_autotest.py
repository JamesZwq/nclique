#!/usr/bin/env python3
# Auto-generated for 5466882

STUDENT_ID = "5466882"
STUDENT_NAME = "Shukun Chen"

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

        n = G.vertex_num
        adj = G.adj_list  # global adjacency list
        tau = [0] * n   # τ_k(v) for every vertex


        for v in range(n):

            Nv = adj[v]
            d  = len(Nv)
            if d == 0:
                continue

            Nv_set = set(Nv)

            loc_adj = {}
            loc_deg = {}
            for u in Nv:
                neigh_u = [w for w in adj[u] if w in Nv_set]
                loc_adj[u] = neigh_u
                loc_deg[u] = len(neigh_u)

            alive = set(Nv)
            peel  = deque([u for u in Nv if loc_deg[u] < k])
            while peel:
                u = peel.popleft()
                if u not in alive:
                    continue
                alive.remove(u)
                for w in loc_adj[u]:
                    if w in alive:
                        loc_deg[w] -= 1
                        if loc_deg[w] < k:
                            peel.append(w)

            comp_cnt = 0
            visited  = set()
            for u in alive:
                if u in visited:
                    continue
                comp_cnt += 1
                stack = [u]
                visited.add(u)
                while stack:
                    x = stack.pop()
                    for y in loc_adj[x]:
                        if y in alive and y not in visited:
                            visited.add(y)
                            stack.append(y)

            tau[v] = comp_cnt      # maybe 0

        return tau


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
