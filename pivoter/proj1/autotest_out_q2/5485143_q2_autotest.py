#!/usr/bin/env python3
# Auto-generated for 5485143

STUDENT_ID = "5485143"
STUDENT_NAME = "Weizhe Yang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # ============== Public API ==============
    @staticmethod
    def process(G, k):
        """
        G : UndirectedUnweightedGraph (vertex ids 0..V-1)
        k : int
        return: List[int]  τ_k(v) for all v
        """
        AL = kCoreBaseStructuralDiversity._extract_adj_list(G)   # list[list[int]], 无向对称
        n = len(AL)
        sd = [0] * n
        for v in range(n):
            sd[v] = kCoreBaseStructuralDiversity._tau_of_vertex(AL, v, k)
        return sd

    # ============== Core ==============
    @staticmethod
    def _tau_of_vertex(AL, v, k):
        N_list = AL[v]
        d = len(N_list)
        if d == 0:
            return 0
        N_set = set(N_list)

        # 为诱导子图 H 构建局部索引
        idx = {u: i for i, u in enumerate(N_list)}
        adjH = [[] for _ in range(d)]
        deg  = [0] * d

        # 建 H 的邻接与度（无向：两端都会计度）
        for u in N_list:
            ui = idx[u]
            for w in AL[u]:
                if w in N_set and w != u:
                    wi = idx[w]
                    adjH[ui].append(wi)
                    deg[ui] += 1

        # k-core 剥皮
        alive = [True] * d
        q = deque([i for i in range(d) if deg[i] < k])
        for i in q:
            alive[i] = False
        while q:
            u = q.popleft()
            for w in adjH[u]:
                if alive[w]:
                    deg[w] -= 1
                    if deg[w] < k:
                        alive[w] = False
                        q.append(w)

        # 统计连通分量
        seen = [False] * d
        comp = 0
        for i in range(d):
            if alive[i] and not seen[i]:
                comp += 1
                dq = deque([i])
                seen[i] = True
                while dq:
                    x = dq.popleft()
                    for y in adjH[x]:
                        if alive[y] and not seen[y]:
                            seen[y] = True
                            dq.append(y)
        return comp

    # ============== Graph extraction ==============
    @staticmethod
    def _extract_adj_list(G):
        """
        输出一个 list[list[int]] 的无向邻接表（保证对称）。
        支持：
          - G.adj_list : list 或 dict
          - G.edges    : list/iter of (u,v) 或 dict项
          - 其他常见字段：adj, neighbors
        """
        # 1) 先尝试 adj_list
        if hasattr(G, "adj_list"):
            return kCoreBaseStructuralDiversity._normalize_adj_list(getattr(G, "adj_list"), G)

        # 2) 再尝试 edges
        if hasattr(G, "edges"):
            return kCoreBaseStructuralDiversity._build_from_edges(getattr(G, "edges"), G)

        # 3) 常见别名
        for name in ("adj", "neighbors", "adjacency", "adj_list_out"):
            if hasattr(G, name):
                return kCoreBaseStructuralDiversity._normalize_adj_list(getattr(G, name), G)

        # 4) 实在没有，尝试空图（避免崩溃）
        n = getattr(G, "vertex_num", 0)
        return [[] for _ in range(n)]

    @staticmethod
    def _normalize_adj_list(adj, G):
        """
        将 adj 正规化为 list[list[int]]，并确保无向对称。
        adj 可以是：
          - list:  每项是邻居列表（元素可能是 int / [v] / (v,...) / {'to':v}）
          - dict:  键为顶点，值为邻居列表
        """
        def parse_vid(x):
            if isinstance(x, int):
                return x
            if isinstance(x, (list, tuple)) and len(x) >= 1 and isinstance(x[0], int):
                return x[0]
            if isinstance(x, dict):
                for key in ("to", "v", "neighbor", "dst", "id"):
                    if key in x and isinstance(x[key], int):
                        return x[key]
            return None

        if isinstance(adj, dict):
            n = getattr(G, "vertex_num", None)
            if n is None:
                if adj:
                    n = max(int(k) for k in adj.keys()) + 1
                else:
                    n = 0
            AL = [[] for _ in range(n)]
            for u, lst in adj.items():
                u = int(u)
                if lst is None: 
                    continue
                for item in lst:
                    v = parse_vid(item)
                    if v is not None:
                        AL[u].append(v)
            # 保证无向对称
            kCoreBaseStructuralDiversity._make_undirected_symmetric(AL)
            return AL

        # list 情况
        if isinstance(adj, (list, tuple)):
            n = len(adj)
            AL = [[] for _ in range(n)]
            for u, lst in enumerate(adj):
                if lst is None:
                    continue
                for item in lst:
                    v = parse_vid(item)
                    if v is not None:
                        AL[u].append(v)
            kCoreBaseStructuralDiversity._make_undirected_symmetric(AL)
            return AL

        # 其他：回退
        n = getattr(G, "vertex_num", 0)
        return [[] for _ in range(n)]

    @staticmethod
    def _build_from_edges(edges, G):
        """
        从边集合构建无向对称邻接表。
        边元素可为 (u,v) / [u,v] / {'u':u,'v':v} / {'src':u,'dst':v} 等。
        """
        def parse_uv(e):
            if isinstance(e, (list, tuple)) and len(e) >= 2:
                a, b = e[0], e[1]
                if isinstance(a, int) and isinstance(b, int):
                    return a, b
            if isinstance(e, dict):
                for ukey, vkey in (("u","v"), ("src","dst"), ("from","to")):
                    if ukey in e and vkey in e and isinstance(e[ukey], int) and isinstance(e[vkey], int):
                        return e[ukey], e[vkey]
            return None

        n = getattr(G, "vertex_num", None)
        max_id = -1
        pairs = []
        for e in edges:
            uv = parse_uv(e)
            if uv is not None:
                u, v = uv
                pairs.append((u, v))
                max_id = max(max_id, u, v)
        if n is None:
            n = max_id + 1 if max_id >= 0 else 0
        AL = [[] for _ in range(n)]
        for u, v in pairs:
            if 0 <= u < n and 0 <= v < n:
                AL[u].append(v)
                AL[v].append(u)
        return AL

    @staticmethod
    def _make_undirected_symmetric(AL):
        """确保无向：如果 u 有 v，则 v 也有 u。"""
        n = len(AL)
        for u in range(n):
            for v in AL[u]:
                if 0 <= v < n and u not in AL[v]:
                    AL[v].append(u)

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
