#!/usr/bin/env python3
# Auto-generated for 5281376

STUDENT_ID = "5281376"
STUDENT_NAME = "Zhongyu Yi"

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
        adj = G.adj_list

        for v in range(n):
            nbrs = adj[v]
            m = len(nbrs)
            if m == 0:
                # no neighbours → zero components
                sd[v] = 0
                continue

            # map each neighbour to an index 0..m-1
            idx = {u: i for i, u in enumerate(nbrs)}
            # build the induced subgraph H on these m nodes
            H_adj = [[] for _ in range(m)]
            for i, u in enumerate(nbrs):
                for w in adj[u]:
                    j = idx.get(w)
                    if j is not None:
                        H_adj[i].append(j)

            # k-core peel: iteratively remove nodes of degree < k
            deg = [len(nei) for nei in H_adj]
            removed = [False] * m
            queue = deque(i for i in range(m) if deg[i] < k)
            rem_count = 0

            while queue:
                u_i = queue.popleft()
                if removed[u_i]:
                    continue
                removed[u_i] = True
                rem_count += 1
                for w_i in H_adj[u_i]:
                    if not removed[w_i]:
                        deg[w_i] -= 1
                        if deg[w_i] < k:
                            queue.append(w_i)

            # if everything was removed, no k-core components
            if rem_count == m:
                sd[v] = 0
                continue

            # count the remaining connected components in H
            visited = [False] * m
            comp_count = 0
            for i in range(m):
                if not removed[i] and not visited[i]:
                    comp_count += 1
                    dq = deque([i])
                    visited[i] = True
                    while dq:
                        x = dq.popleft()
                        for y in H_adj[x]:
                            if not removed[y] and not visited[y]:
                                visited[y] = True
                                dq.append(y)

            sd[v] = comp_count

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
