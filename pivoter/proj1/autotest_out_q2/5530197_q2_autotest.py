#!/usr/bin/env python3
# Auto-generated for 5530197

STUDENT_ID = "5530197"
STUDENT_NAME = "Fengyuan Liu"

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
        adj = G.adj_list
        tau = [0] * n

        if k < 0:  # Illegal k, treat it as 0
            return tau

        # mark[v] == tag means v belongs to the neighbor set of the vertex currently being processed
        mark = [-1] * n
        tag = 0

        for v in range(n):
            nbrs = adj[v]
            d_v = len(nbrs)

            # Empty neighbors or no k-core possible (maximum degree ≤ d_v-1)
            if d_v == 0:
                tau[v] = 0
                continue
            if k > d_v - 1:
                tau[v] = 0
                continue

            # Mark neighbor set
            tag += 1
            local_index = {}
            for i, u in enumerate(nbrs):
                mark[u] = tag
                local_index[u] = i

            # Construct neighbor-induced subgraph local adjacency list (0..d_v-1)
            local_adj = [[] for _ in range(d_v)]
            for i_u, u in enumerate(nbrs):
                for w in adj[u]:
                    if mark[w] == tag:
                        i_w = local_index[w]
                        if i_u < i_w:            # Avoid duplicate insertion
                            local_adj[i_u].append(i_w)
                            local_adj[i_w].append(i_u)

            # k-core peeling
            keep = kCoreBaseStructuralDiversity._kcore_peel(local_adj, k)

            if not keep or not any(keep):
                tau[v] = 0
            else:
                tau[v] = kCoreBaseStructuralDiversity._count_components(local_adj, keep)

        return tau
    
    @staticmethod
    def _kcore_peel(local_adj, k):
        """
        Queue-based k-core peeling on a *local* undirected adjacency list.
        """
        m = len(local_adj)
        if m == 0:
            return []
        if k <= 0:
            return [True] * m

        deg = [len(local_adj[i]) for i in range(m)]
        keep = [True] * m
        q = deque([i for i, d in enumerate(deg) if d < k])

        while q:
            u = q.popleft()
            if not keep[u]:
                continue
            keep[u] = False
            for w in local_adj[u]:
                if keep[w]:
                    deg[w] -= 1
                    if deg[w] == k - 1:  # just dropped below k
                        q.append(w)
        return keep

    @staticmethod
    def _count_components(local_adj, keep):
        """
        Count connected components of the subgraph induced by nodes with keep[i]=True.
        """
        m = len(local_adj)
        visited = [False] * m
        comps = 0
        for i in range(m):
            if not keep[i] or visited[i]:
                continue
            comps += 1
            dq = deque([i])
            visited[i] = True
            while dq:
                u = dq.popleft()
                for w in local_adj[u]:
                    if keep[w] and not visited[w]:
                        visited[w] = True
                        dq.append(w)
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
