#!/usr/bin/env python3
# Auto-generated for 5451986

STUDENT_ID = "5451986"
STUDENT_NAME = "Xuewei Yao"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        # Time Complexity: O(n * d_max^2)
        n = G.vertex_num
        tau = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]

            if len(neighbors) < k:
                tau[v] = 0
                continue

            # Build neighbor-induced subgraph G[N(v)]
            neighbor_set = set(neighbors)
            neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbors)}

            induced_adj = [[] for _ in range(len(neighbors))]
            for i, u in enumerate(neighbors):
                for w in G.adj_list[u]:
                    if w in neighbor_set:
                        induced_adj[i].append(neighbor_to_idx[w])

            tau[v] = kCoreBaseStructuralDiversity.count_k_cores(induced_adj, k)

        return tau

    @staticmethod
    def count_k_cores(adj_list, k):
        n = len(adj_list)
        if n < k:
            return 0

        degrees = [len(adj_list[i]) for i in range(n)]
        active = [True] * n

        # Remove vertices with degree < k
        queue = deque()
        for v in range(n):
            if degrees[v] < k:
                queue.append(v)

        while queue:
            v = queue.popleft()
            if not active[v]:
                continue
            active[v] = False

            for u in adj_list[v]:
                if active[u]:
                    degrees[u] -= 1
                    if degrees[u] < k:
                        queue.append(u)

        # Count connected components
        visited = [False] * n
        num_components = 0

        for v in range(n):
            if active[v] and not visited[v]:
                component_queue = deque([v])
                visited[v] = True

                while component_queue:
                    u = component_queue.popleft()
                    for w in adj_list[u]:
                        if active[w] and not visited[w]:
                            visited[w] = True
                            component_queue.append(w)

                num_components += 1

        return num_components

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
