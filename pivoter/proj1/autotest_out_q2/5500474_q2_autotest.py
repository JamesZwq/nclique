#!/usr/bin/env python3
# Auto-generated for 5500474

STUDENT_ID = "5500474"
STUDENT_NAME = "Jiayi Sun"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(graph, k):
        """
        Compute τ_k(v) for all nodes v in the graph.
        τ_k(v) is the number of connected components in the k-core
        of the subgraph induced by neighbors of v.

        Parameters
        ----------
        graph : UndirectedUnweightedGraph
        k : int

        Returns
        -------
        List[int] : τ_k(v) for all v

        Time Complexity:
        ----------------
        Let d_v be the degree of node v.
        For each node:
            - Neighbor collection: O(d_v)
            - Subgraph construction: O(d_v^2)
            - k-core pruning: O(d_v^2)
            - DFS for connected components: O(d_v^2)
        Worst-case total: O(n^3) for dense graphs.
        Best-case total: O(n) for sparse graphs.
        """
        n = graph.vertex_num
        res = [0] * n

        for v in range(n):
            neigh = set(graph.adj_list[v])
            if not neigh:
                continue

            # Subgraph construction (induced by neighbors of v)
            # Time: O(d_v^2)
            sub_g = {u: [] for u in neigh}
            for u in neigh:
                sub_g[u] = [w for w in graph.adj_list[u] if w in neigh]

            # Count k-core connected components in the subgraph
            res[v] = kCoreBaseStructuralDiversity._count_kcore_components(sub_g, k)

        return res

    @staticmethod
    def _count_kcore_components(sub_g, k):
        """
        Compute number of connected components in the k-core of sub_g.

        Parameters
        ----------
        sub_g : Dict[int, List[int]]
            An undirected subgraph represented as adjacency list.
        k : int
            Minimum degree for k-core.

        Returns
        -------
        int : Number of connected components in the k-core.

        Time Complexity:
        ----------------
        - Degree initialization: O(d^2)
        - k-core pruning: O(d^2)
        - DFS traversal: O(d^2)
        where d is the number of nodes in sub_g (i.e., d <= d_v).
        """
        if not sub_g:
            return 0

        # Step 1: Initialize degree of each node — O(d)
        deg = {u: len(sub_g[u]) for u in sub_g}
        drop = set()
        changed = True

        # Step 2: Iteratively remove nodes with degree < k — O(d^2) worst-case
        while changed:
            changed = False
            for u in sub_g:
                if u in drop:
                    continue
                if deg[u] < k:
                    drop.add(u)
                    changed = True
                    for v in sub_g[u]:
                        if v not in drop:
                            deg[v] -= 1

        # Step 3: Remaining nodes form k-core — O(d)
        kept = [u for u in sub_g if u not in drop]
        if not kept:
            return 0

        # Step 4: Count connected components via DFS — O(d^2) in dense case
        seen = set()
        cnt = 0

        for start in kept:
            if start in seen:
                continue

            cnt += 1
            stk = [start]
            idx = 0

            while idx < len(stk):
                cur = stk[idx]
                idx += 1
                if cur in seen:
                    continue
                seen.add(cur)
                push_n = [w for w in sub_g[cur] if w in kept and w not in seen]
                if push_n:
                    stk.extend(push_n)

        return cnt

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
