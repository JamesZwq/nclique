#!/usr/bin/env python3
# Auto-generated for 5356286

STUDENT_ID = "5356286"
STUDENT_NAME = "Yifan Yang"

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
        n = G.vertex_num

        degree = kCoreBaseStructuralDiversity.getDegree(G)

        degeneracyOrder, coreNum = \
            kCoreBaseStructuralDiversity.coreDecomposition(G, degree)

        newG = \
            kCoreBaseStructuralDiversity.getOrientedGraph(G, degeneracyOrder)

        bitmap = [0] * n
        stack = deque()

        sd = [0] * n

        for v in range(n):
            inducedGraph = [list() for _ in range(n)]

            for u in G.adj_list[v]:
                bitmap[u] = 1
                stack.append(u)

            for u in G.adj_list[v]:
                for w in newG[u]:
                    if (coreNum[w] >= k):
                        if (bitmap[w] == 1):
                            inducedGraph[u].append(w)
                            inducedGraph[w].append(u)

            sd[v] = kCoreBaseStructuralDiversity.kCore(inducedGraph, k, bitmap)

            while (stack):
                u = stack.pop()
                bitmap[u] = 0

        return sd


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def getDegree(G):
        n = G.vertex_num
        ans = [0] * n

        for u in range(n):
            ans[u] = len(G.adj_list[u])
        return ans

    @staticmethod
    def coreDecomposition(G, degree):
        n = G.vertex_num

        d = degree[:]
        b = [0] * n
        D, p = kCoreBaseStructuralDiversity.binSort(degree, b)

        degeneracyOrder = [0] * n

        for i in range(n):
            v = D[i]
            degeneracyOrder[v] = i
            for u in G.adj_list[v]:
                if (d[u] > d[v]):
                    du = d[u]
                    pu = p[u]
                    pw = b[du]
                    w = D[pw]
                    if (u != w):
                        D[pu] = w
                        D[pw] = u
                        p[u] = pw
                        p[w] = pu
                    b[du] = b[du] + 1
                    d[u] = d[u] - 1

        return (degeneracyOrder, d)

    @staticmethod
    def binSort(degree, b):
        n = len(degree)
        bins = [[] for _ in range(n)]

        for i in range(n):
            bins[degree[i]].append((i))

        sortedList = [0] * n
        position = [0] * n

        ansCnt = 0
        bcnt = 0
        for bin in bins:
            b[bcnt] = ansCnt
            bcnt = bcnt + 1

            for i in range(len(bin)):
                sortedList[ansCnt] = bin[i]
                position[bin[i]] = ansCnt
                ansCnt = ansCnt + 1
        return (sortedList, position)

    @staticmethod
    def getOrientedGraph(G, degeneracyOrder):
        n = G.vertex_num
        newG = [list() for _ in range(n)]

        for u in range(n):
            for v in G.adj_list[u]:
                if (degeneracyOrder[u] < degeneracyOrder[v]):
                    newG[u].append(v)
        return newG

    @staticmethod
    def kCore(G, k, bitmap):
        n = len(G)
        degree = [0] * n
        visited = [False] * n

        for u in range(n):
            if bitmap[u]:
                degree[u] = sum(1 for v in G[u] if bitmap[v])

        queue = deque()
        for i in range(n):
            if bitmap[i] and degree[i] < k:
                queue.append(i)
                visited[i] = True

        while queue:
            u = queue.popleft()
            for v in G[u]:
                if bitmap[v] and not visited[v]:
                    degree[v] -= 1
                    if degree[v] < k:
                        visited[v] = True
                        queue.append(v)

        def bfs(start):
            q = deque([start])
            visited[start] = True
            while q:
                u = q.popleft()
                for v in G[u]:
                    if bitmap[v] and not visited[v] and degree[v] >= k:
                        visited[v] = True
                        q.append(v)

        count = 0
        for i in range(n):
            if bitmap[i] and not visited[i] and degree[i] >= k:
                bfs(i)
                count += 1

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
