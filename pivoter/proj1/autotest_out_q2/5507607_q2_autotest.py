#!/usr/bin/env python3
# Auto-generated for 5507607

STUDENT_ID = "5507607"
STUDENT_NAME = "Yuxi Chang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # Compute the global k-core membership of each vertex
    @staticmethod
    def global_k_core(G, k):
        n = G.vertex_num
        degrees = [len(G.adj_list[v]) for v in range(n)]
        queue = deque([v for v in range(n) if degrees[v] < k])
        in_core = [True] * n

        while queue:
            u = queue.popleft()
            in_core[u] = False
            for w in G.adj_list[u]:
                if in_core[w]:
                    degrees[w] -= 1
                    if degrees[w] < k:
                        queue.append(w)
        return in_core

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
        sd = [0] * n

        # Compute global k-core membership for all vertices
        in_core = kCoreBaseStructuralDiversity.global_k_core(G, k)

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                continue

            # Filter neighbors to keep only those in the global k-core
            neighbor_set = {u for u in neighbors if in_core[u]}
            if not neighbor_set:
                continue

            # Build the induced subgraph for these neighbors
            sub_adj = {u: [w for w in G.adj_list[u] if w in neighbor_set] for u in neighbor_set}

            # Perform local peeling: remove nodes with degree less than k
            degrees = {u: len(sub_adj[u]) for u in neighbor_set}
            queue = deque([u for u in degrees if degrees[u] < k])
            while queue:
                u = queue.popleft()
                for w in sub_adj[u]:
                    if w in degrees:
                        degrees[w] -= 1
                        if degrees[w] == k - 1:
                            queue.append(w)
                del degrees[u]

            # Count connected components in the remaining subgraph
            visited = set()
            count = 0
            for u in degrees:
                if u not in visited:
                    q = deque([u])
                    visited.add(u)
                    while q:
                        x = q.popleft()
                        for y in sub_adj[x]:
                            if y in degrees and y not in visited:
                                visited.add(y)
                                q.append(y)
                    count += 1

            sd[v] = count

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
