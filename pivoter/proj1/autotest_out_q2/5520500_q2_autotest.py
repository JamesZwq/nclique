#!/usr/bin/env python3
# Auto-generated for 5520500

STUDENT_ID = "5520500"
STUDENT_NAME = "Yixuan Liu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def map_neighbors_to_indices(neighbors):
        idx = {}
        for i, u in enumerate(neighbors):
            idx[u] = i
        return idx

    @staticmethod
    def build_local_adjacency(G, neighbors, idx):
        m = len(neighbors)
        local_adj = [[] for _ in range(m)]
        for u in neighbors:
            iu = idx[u]
            for w in G.adj_list[u]:
                if w in idx:
                    local_adj[iu].append(idx[w])
        return local_adj

    @staticmethod
    def peeling_k_core(local_adj, deg_local, k):
        m = len(local_adj)
        removed = [False] * m
        q = deque(i for i in range(m) if deg_local[i] < k)

        while q:
            u = q.popleft()
            if removed[u]:
                continue
            removed[u] = True
            for w in local_adj[u]:
                if not removed[w]:
                    deg_local[w] -= 1
                    if deg_local[w] == k - 1:
                        q.append(w)
        return removed

    @staticmethod
    def count_components(local_adj, removed):
        m = len(local_adj)
        visited = [False] * m
        comp_count = 0

        for i in range(m):
            if not removed[i] and not visited[i]:
                comp_count += 1
                bfs_queue = deque([i])
                visited[i] = True
                while bfs_queue:
                    x = bfs_queue.popleft()
                    for y in local_adj[x]:
                        if not removed[y] and not visited[y]:
                            visited[y] = True
                            bfs_queue.append(y)
        return comp_count

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            m = len(neighbors)
            if m == 0:
                sd[v] = 0
                continue

            idx = kCoreBaseStructuralDiversity.map_neighbors_to_indices(neighbors)
            local_adj = kCoreBaseStructuralDiversity.build_local_adjacency(G, neighbors, idx)
            deg_local = [len(neis) for neis in local_adj]
            removed = kCoreBaseStructuralDiversity.peeling_k_core(local_adj, deg_local, k)
            comp_count = kCoreBaseStructuralDiversity.count_components(local_adj, removed)
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
