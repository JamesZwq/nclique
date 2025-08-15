#!/usr/bin/env python3
# Auto-generated for 5461155

STUDENT_ID = "5461155"
STUDENT_NAME = "Yifei Li"

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
        adj = G.adj_list
        flag = [-1] * n  # Mark neighbors
        pruned = [-1] * n  # Track removed vertices
        visited = [-1] * n  # Track visited vertices
        epoch = 0  # Timestamp for marking

        tau = [0] * n

        for v in range(n):
            nbrs = adj[v]
            deg_v = len(nbrs)

            if deg_v < k + 1:  # Skip if not enough neighbors
                continue

            epoch += 1
            kCoreBaseStructuralDiversity._flag_neighbors(nbrs, flag, epoch)

            if k >= 1:
                deg = kCoreBaseStructuralDiversity._calc_induced_degrees(adj, nbrs, flag, epoch)
                kCoreBaseStructuralDiversity._peel_kcore(adj, nbrs, deg, k, flag, pruned, epoch)

            comps = kCoreBaseStructuralDiversity._count_components(adj, nbrs, flag, pruned, visited, epoch)
            tau[v] = comps

        return tau

    @staticmethod
    def _flag_neighbors(nbrs, flag, epoch):
        """Mark neighbors with current epoch."""
        for u in nbrs:
            flag[u] = epoch

    @staticmethod
    def _calc_induced_degrees(adj, nbrs, flag, epoch):
        """Calculate degrees in the induced subgraph."""
        deg = {}
        for u in nbrs:
            cnt = sum(1 for w in adj[u] if flag[w] == epoch)
            deg[u] = cnt
        return deg

    @staticmethod
    def _peel_kcore(adj, nbrs, deg, k, flag, pruned, epoch):
        """Peel vertices not meeting k-core criteria."""
        to_check = [u for u in nbrs if deg[u] < k]
        head = 0

        while head < len(to_check):
            u = to_check[head]
            head += 1

            if pruned[u] == epoch:
                continue
            pruned[u] = epoch

            neighbours = adj[u]
            for w in neighbours:
                if flag[w] != epoch or pruned[w] == epoch:
                    continue
                deg[w] -= 1
                if deg[w] == k - 1:
                    to_check.append(w)

    @staticmethod
    def _count_components(adj, nbrs, flag, pruned, visited, epoch):
        """Count connected components in the remaining graph."""
        comps = 0
        for src in nbrs:
            if pruned[src] == epoch or visited[src] == epoch:
                continue

            comps += 1
            q = deque([src])
            visited[src] = epoch
            while len(q) > 0:
                u = q.popleft()
                neighbours = adj[u]
                for w in neighbours:
                    if flag[w] != epoch or pruned[w] == epoch or visited[w] == epoch:
                        continue
                    visited[w] = epoch
                    q.append(w)

        return comps

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
