#!/usr/bin/env python3
# Auto-generated for 5467882

STUDENT_ID = "5467882"
STUDENT_NAME = "Xi He"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k, strategy='minscan'):
        """
        返回数组 sd，其中 sd[v] = τ_k(v)
        strategy:
            'naive'    -> 全部使用朴素构建
            'minscan'  -> 全部使用 min-scan
            'auto'     -> 顶点级别自动选择
        """
        # 统一邻接表
        if hasattr(G, 'adj_list'):
            adj = G.adj_list
        elif hasattr(G, 'adj_list_out'):
            adj = G.adj_list_out
        else:
            raise AttributeError("Graph must have adj_list or adj_list_out")

        n = G.vertex_num
        sd = [0] * n

        # 在需要 min-scan 时才构建邻接集合
        need_adj_set = (strategy != 'naive')
        adj_set = [set(lst) for lst in adj] if need_adj_set else None

        # 复用的辅助结构
        mark = [0] * n
        current_tag = 1
        deg_in_sub = [0] * n
        removed = [False] * n
        visited = [False] * n
        sub_neighbors = [[] for _ in range(n)]

        for v in range(n):
            neighbors = adj[v]
            d_v = len(neighbors)
            if d_v == 0:
                sd[v] = 0
                continue

            current_tag += 1
            # 初始化邻居状态
            for u in neighbors:
                mark[u] = current_tag
                deg_in_sub[u] = 0
                removed[u] = False
                visited[u] = False
                sub_neighbors[u].clear()

            # ---------- 选择策略 ----------
            local_strategy = strategy
            if strategy == 'auto':
                # 估算两种方法代价
                sum_deg = 0
                sum_min = 0
                for u in neighbors:
                    du = len(adj[u])
                    sum_deg += du
                    sum_min += du if du <= d_v else d_v
                # 如果朴素更便宜或一样，则用朴素
                local_strategy = 'naive' if sum_deg <= sum_min else 'minscan'

            # ---------- 构建邻居诱导子图 ----------
            if local_strategy == 'naive':
                # 朴素：对每个邻居扫描其全部邻接表
                for u in neighbors:
                    for w in adj[u]:
                        if w == v or mark[w] != current_tag:
                            continue
                        if u < w:  # 无向边去重
                            sub_neighbors[u].append(w)
                            sub_neighbors[w].append(u)
                            deg_in_sub[u] += 1
                            deg_in_sub[w] += 1
            else:
                # min-scan
                for u in neighbors:
                    du = len(adj[u])
                    if du <= d_v:
                        for w in adj[u]:
                            if w == v or mark[w] != current_tag:
                                continue
                            if u < w:
                                sub_neighbors[u].append(w)
                                sub_neighbors[w].append(u)
                                deg_in_sub[u] += 1
                                deg_in_sub[w] += 1
                    else:
                        uset = adj_set[u]
                        for w in neighbors:
                            if w == u:
                                continue
                            if u < w and w in uset:
                                sub_neighbors[u].append(w)
                                sub_neighbors[w].append(u)
                                deg_in_sub[u] += 1
                                deg_in_sub[w] += 1

            # ---------- k-core 剥除 ----------
            if k > 0:
                q = deque()
                for u in neighbors:
                    if deg_in_sub[u] < k:
                        removed[u] = True
                        q.append(u)
                while q:
                    x = q.popleft()
                    for w in sub_neighbors[x]:
                        if not removed[w]:
                            deg_in_sub[w] -= 1
                            if deg_in_sub[w] < k:
                                removed[w] = True
                                q.append(w)

            # ---------- 统计连通分量 ----------
            comp = 0
            for u in neighbors:
                if removed[u] or visited[u]:
                    continue
                comp += 1
                dq = deque([u])
                visited[u] = True
                while dq:
                    x = dq.popleft()
                    for w in sub_neighbors[x]:
                        if not removed[w] and not visited[w]:
                            visited[w] = True
                            dq.append(w)
            sd[v] = comp

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
