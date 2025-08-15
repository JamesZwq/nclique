#!/usr/bin/env python3
# Auto-generated for 5538454

STUDENT_ID = "5538454"
STUDENT_NAME = "Yuting Yao"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        for i in range(n):
            pals = []
            for m in G.adj_list[i]:
                pals.append(m)
            if len(pals) == 0:
                sd[i] = 0
                continue
            ref = {}
            for p in range(len(pals)):
                ref[pals[p]] = p
            N = len(pals)
            adj = []
            for _ in range(N):
                adj.append([])
            for u in range(N):
                node = pals[u]
                links = G.adj_list[node]
                for to in links:
                    if (to in ref) and (to != node):
                        adj[u].append(ref[to])
            alive = []
            for _ in range(N):
                alive.append(1)
            degrees = []
            for x in range(N):
                degrees.append(len(adj[x]))

            found = 1
            while found == 1:
                found = 0
                for p in range(N):
                    if alive[p] == 1 and degrees[p] < k:
                        alive[p] = 0
                        found = 1
                        for q in adj[p]:
                            if alive[q] == 1:
                                degrees[q] = degrees[q] - 1
            tagged = []
            for _ in range(N):
                tagged.append(0)
            num = 0
            for start in range(N):
                if alive[start] == 1 and tagged[start] == 0:
                    que = [start]
                    tagged[start] = 1
                    ptr = 0
                    while ptr < len(que):
                        take = que[ptr]
                        ptr = ptr + 1
                        edges = adj[take]
                        for point in edges:
                            if alive[point] == 1 and tagged[point] == 0:
                                tagged[point] = 1
                                que.append(point)
                    num = num + 1
            sd[i] = num
        return sd
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
