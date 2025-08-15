#!/usr/bin/env python3
# Auto-generated for 5511254

STUDENT_ID = "5511254"
STUDENT_NAME = "David Cai"

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
        # number of vertices in G
        total_verts = G.vertex_num                      # O(1)
        # initialize result array
        diversity = [0] * total_verts                   # O(n)
        # adjacency list reference
        adj_map = G.adj_list                             # O(1)

        # iterate over each vertex
        v = 0                                            # O(1)
        while v < total_verts:                          # O(n) iterations
            # build neighbor set of v
            nbrs_v = set(adj_map[v])                    # O(d_v)
            # build induced subgraph H on neighbors
            H = {}                                      # O(1)
            for u in nbrs_v:                            # O(d_v)
                # intersect u's neighbors with nbrs_v
                H[u] = set(adj_map[u]) & nbrs_v         # O(d_u + d_v)

            # initialize peel queue with nodes of degree < k
            queue = deque()                             # O(1)
            for u, nbrs in H.items():                   # O(|nbrs_v|)
                if len(nbrs) < k:                       # O(1)
                    queue.append(u)                     # O(1)

            # peel low‑degree nodes
            while queue:                                # O(|nbrs_v|) total pops
                u0 = queue.popleft()                    # O(1)
                rem = H.pop(u0, None)                   # O(1)
                if rem is None:                         # O(1)
                    continue                            # O(1)
                for w in rem:                           # O(deg_H(u0))
                    if w in H:                          # O(1)
                        H[w].remove(u0)                 # O(1)
                        if len(H[w]) < k:               # O(1)
                            queue.append(w)             # O(1)

            # count connected components in remaining H
            seen = set()                                # O(1)
            comps = 0                                   # O(1)
            for start in H:                             # O(|H|)
                if start not in seen:                  # O(1)
                    comps += 1                         # O(1)
                    bfs = deque([start])               # O(1)
                    seen.add(start)                    # O(1)
                    while bfs:                         # O(|H| + edges in H)
                        cur = bfs.popleft()            # O(1)
                        for w in H[cur]:               # O(deg_H(cur))
                            if w not in seen:         # O(1)
                                seen.add(w)            # O(1)
                                bfs.append(w)          # O(1)

            # record the component count for v
            diversity[v] = comps                        # O(1)
            v += 1                                      # O(1)

        return diversity                                # O(1)
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
