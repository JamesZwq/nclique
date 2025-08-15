#!/usr/bin/env python3
# Auto-generated for 5554005

STUDENT_ID = "5554005"
STUDENT_NAME = "Xiaoyang Chen"

# ======= 学生代码 =======
from collections import deque

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
        τ = [0] * G.vertex_num
        for v in range(G.vertex_num):
            # 1) 构建邻居诱导子图的动态邻接表 subG
            nbr = G.adj_list[v]
            if not nbr:
                continue
            subG = {u: set() for u in nbr}
            for u in nbr:
                for w in G.adj_list[u]:
                    if w in subG:
                        subG[u].add(w)
            # 2) k-core 修剪：不断删除度 < k 的点
            q = deque(u for u, neigh in subG.items() if len(neigh) < k)
            while q:
                u0 = q.popleft()
                if u0 not in subG:
                    continue
                for w in subG[u0]:
                    subG[w].remove(u0)
                    if len(subG[w]) < k:
                        q.append(w)
                del subG[u0]

            # 3) 剩余 subG 中跑一次连通分量
            visited = set()
            cnt = 0
            for u in subG:
                if u in visited:
                    continue
                cnt += 1
                stack = [u]
                visited.add(u)

                while stack:
                    x = stack.pop()
                    for w in subG[x]:
                        if w not in visited:
                            visited.add(w)
                            stack.append(w)

            τ[v] = cnt

        return τ

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
