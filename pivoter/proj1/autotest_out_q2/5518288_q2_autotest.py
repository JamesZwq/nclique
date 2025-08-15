#!/usr/bin/env python3
# Auto-generated for 5518288

STUDENT_ID = "5518288"
STUDENT_NAME = "Bingyu Wang"

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
        if n == 0:
            return []

        # 将邻接列表转换为集合以提高查询效率
        adj_set = [set(neighbors) for neighbors in G.adj_list]
        tau = [0] * n

        for v in range(n):
            neighbors_v = G.adj_list[v]
            d_v = len(neighbors_v)

            # 如果顶点没有邻居，则τ_k(v)=0
            if d_v == 0:
                continue

            S = set(neighbors_v)  # 邻居集合
            sub_adj = {}  # 子图的邻接表

            # 构建邻居诱导子图
            for u in neighbors_v:
                # 获取u的邻居与S的交集
                common = adj_set[u] & S
                sub_adj[u] = list(common)

            # 初始化度数和删除标记
            degree = {u: len(sub_adj[u]) for u in S}
            deleted = {u: False for u in S}
            q = deque()

            # 初始化队列：添加所有度小于k的顶点
            for u in S:
                if degree[u] < k:
                    deleted[u] = True
                    q.append(u)

            # k-core分解过程
            while q:
                u = q.popleft()
                for w in sub_adj[u]:
                    if not deleted[w]:
                        degree[w] -= 1
                        if degree[w] < k:
                            deleted[w] = True
                            q.append(w)

            # 收集未删除的顶点
            remaining = [u for u in S if not deleted[u]]

            # 如果没有剩余顶点，τ_k(v)=0
            if not remaining:
                continue

            # 计算剩余子图的连通分量数量
            visited = {u: False for u in remaining}
            count = 0

            for u in remaining:
                if not visited[u]:
                    count += 1
                    stack = [u]
                    visited[u] = True
                    while stack:
                        cur = stack.pop()
                        for neighbor in sub_adj[cur]:
                            # 跳过已删除或已访问的邻居
                            if deleted[neighbor] or visited[neighbor]:
                                continue
                            visited[neighbor] = True
                            stack.append(neighbor)

            tau[v] = count

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
