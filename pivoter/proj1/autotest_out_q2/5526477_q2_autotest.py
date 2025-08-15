#!/usr/bin/env python3
# Auto-generated for 5526477

STUDENT_ID = "5526477"
STUDENT_NAME = "(Eric) Qinfeng Luo"

# ======= 学生代码 =======
from collections import deque
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):

        n = G.vertex_num
        result = [0] * n

        # Traverse all vertices to generate neighbor-induced subgraph
        for v in range(n):
            # neighbors of vertex v
            neighbors = G.adj_list[v]

            if len(neighbors) == 0:
                result[v] = 0
                continue

            # Generate neighbor-induced subgraph
            neighbor_set = set(neighbors)
            neighbor_list_removeDuplicate = list(dict.fromkeys(neighbors))

            neighbor_idx = {neighbor: i for i, neighbor in enumerate(neighbor_list_removeDuplicate)}
            adj_list_neighbor = [[] for _ in range(len(neighbor_list_removeDuplicate))]
            for i, u in enumerate(neighbor_list_removeDuplicate):
                for w in G.adj_list[u]:
                    # w>u is to control no duplicate edge is added
                    if w in neighbor_set and w > u:
                        j = neighbor_idx[w]
                        adj_list_neighbor[i].append(j)
                        adj_list_neighbor[j].append(i)

            # Count k-cores in the neighbor-induced subgraph
            result[v] = kCoreBaseStructuralDiversity.k_cores(adj_list_neighbor, k)

        return result

    # Count the number of k-core components in the sub graph
    @staticmethod
    def k_cores(adj_list, k):
        n = len(adj_list)
        if n == 0:
            return 0

        # Find k-core using Peeling
        degrees = [len(adj_list[i]) for i in range(n)]
        removed = [False] * n
        queue = deque()

        # adding vertices into queue when degree < k
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)
                removed[i] = True

        # remove vertices with degree < k
        while queue:
            u = queue.popleft()
            for v in adj_list[u]:
                if not removed[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)
                        removed[v] = True

        # Count connected components using DFS
        visited = [False] * n
        components = 0

        for node in range(n):
            if not removed[node] and not visited[node]:
                kCoreBaseStructuralDiversity.dfs(node, adj_list, visited, removed)
                components += 1

        return components

    #  DFS to explore connected component in k-core
    @staticmethod
    def dfs(start_node, adj_list, visited, removed):
        stack = [start_node]
        visited[start_node] = True

        while stack:
            u = stack.pop()
            for v in adj_list[u]:
                if not removed[v] and not visited[v]:
                    visited[v] = True
                    stack.append(v)

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
