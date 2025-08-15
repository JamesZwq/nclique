#!/usr/bin/env python3
# Auto-generated for 5531741

STUDENT_ID = "5531741"
STUDENT_NAME = "Chengbo Zhang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        adj = G.adj_list
        tau_k = [0] * n  # output array

        for v in range(n):
            neighbors = adj[v]
            if not neighbors:
                continue

            # Step 1: Build neighbor-induced subgraph nbr(v)
            nbr_nodes = list(set(neighbors))  # no self-loop
            nbr_id_map = {node: idx for idx, node in enumerate(nbr_nodes)}  # mapping to local indices
            m = len(nbr_nodes)
            nbr_adj = [[] for _ in range(m)]

            for i, u in enumerate(nbr_nodes):
                for w in adj[u]:
                    if w in nbr_id_map and w != v:  # avoid self-loop or edges to v
                        nbr_adj[i].append(nbr_id_map[w])

            # Step 2: Compute k-core inside neighbor-induced subgraph
            deg = [len(nbr_adj[i]) for i in range(m)]
            removed = [False] * m
            q = deque([i for i in range(m) if deg[i] < k])

            while q:
                u = q.popleft()
                if removed[u]:
                    continue
                removed[u] = True
                for w in nbr_adj[u]:
                    if not removed[w]:
                        deg[w] -= 1
                        if deg[w] == k - 1:
                            q.append(w)

            # Step 3: Count number of connected components in remaining subgraph
            visited = [False] * m
            count = 0

            for i in range(m):
                if removed[i] or visited[i]:
                    continue
                count += 1
                q = deque([i])
                visited[i] = True
                while q:
                    u = q.popleft()
                    for w in nbr_adj[u]:
                        if not removed[w] and not visited[w]:
                            visited[w] = True
                            q.append(w)

            tau_k[v] = count

        return tau_k

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
