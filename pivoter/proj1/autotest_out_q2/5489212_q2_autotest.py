#!/usr/bin/env python3
# Auto-generated for 5489212

STUDENT_ID = "5489212"
STUDENT_NAME = "Bingfeng Lin"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        def get_k_core_subgraph(graph_adj, k):
            n = len(graph_adj)
            degree = [len(neigh) for neigh in graph_adj]
            removed = [False] * n
            queue = deque()

            for i in range(n):
                if degree[i] < k:
                    queue.append(i)
                    removed[i] = True

            while queue:
                v = queue.popleft()
                for u in graph_adj[v]:
                    if not removed[u]:
                        degree[u] -= 1
                        if degree[u] < k:
                            removed[u] = True
                            queue.append(u)

            # build subgraph
            subgraph = [[] for _ in range(n)]
            for u in range(n):
                if not removed[u]:
                    for v in graph_adj[u]:
                        if not removed[v]:
                            subgraph[u].append(v)
            return subgraph, removed

        def count_connected_components(graph_adj, removed):
            n = len(graph_adj)
            visited = [False] * n
            count = 0

            for i in range(n):
                if not removed[i] and not visited[i]:
                    count += 1
                    queue = deque([i])
                    visited[i] = True
                    while queue:
                        u = queue.popleft()
                        for v in graph_adj[u]:
                            if not visited[v]:
                                visited[v] = True
                                queue.append(v)
            return count

        n = G.vertex_num
        tau_k = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                tau_k[v] = 0
                continue

            node_list = list(set(neighbors))
            id_map = {old: new for new, old in enumerate(node_list)}
            subgraph = [[] for _ in range(len(node_list))]

            for u in node_list:
                for w in G.adj_list[u]:
                    if w in id_map and w != v:
                        subgraph[id_map[u]].append(id_map[w])

            kcore_subgraph, removed = get_k_core_subgraph(subgraph, k)
            tau_k[v] = count_connected_components(kcore_subgraph, removed)

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
