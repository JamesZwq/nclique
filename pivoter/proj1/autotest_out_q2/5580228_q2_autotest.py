#!/usr/bin/env python3
# Auto-generated for 5580228

STUDENT_ID = "5580228"
STUDENT_NAME = "(Zaine) Fan Liu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            nbrs = G.adj_list[v]
            m = len(nbrs)
            if m == 0:
                # No neighbors, diversity is zero
                sd[v] = 0
                continue

            # Map neighbor vertex to index in induced subgraph
            idx = {u: i for i, u in enumerate(nbrs)}

            # Build induced adjacency and degree arrays
            deg = [0] * m
            adj = [[] for _ in range(m)]
            for i, u in enumerate(nbrs):
                for w in G.adj_list[u]:
                    j = idx.get(w)
                    if j is not None:
                        deg[i] += 1
                        adj[i].append(j)

            # Peeling: remove vertices with degree < k iteratively
            queue = deque([i for i in range(m) if deg[i] < k])
            removed = [False] * m
            while queue:
                u = queue.popleft()
                if removed[u]:
                    continue
                removed[u] = True
                for w in adj[u]:
                    if not removed[w]:
                        deg[w] -= 1
                        if deg[w] == k - 1:
                            queue.append(w)

            # Count connected components among remaining vertices
            visited = [False] * m
            comps = 0
            for i in range(m):
                if not removed[i] and not visited[i]:
                    comps += 1
                    # DFS/BFS to mark this component
                    stack = [i]
                    visited[i] = True
                    while stack:
                        x = stack.pop()
                        for y in adj[x]:
                            if not removed[y] and not visited[y]:
                                visited[y] = True
                                stack.append(y)

            sd[v] = comps

        return sd

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
