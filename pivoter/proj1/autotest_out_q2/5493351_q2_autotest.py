#!/usr/bin/env python3
# Auto-generated for 5493351

STUDENT_ID = "5493351"
STUDENT_NAME = "Jintao Zhang"

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
        tau = [0] * n

        for v in range(n):
            neighbors = set(G.adj_list[v])

            #handle vertices with no neighbors
            if len(neighbors) == 0:
                tau[v] = 0
                continue

            #compute k-core diversity for this vertex
            tau[v] = kCoreBaseStructuralDiversity._compute_k_core_diversity(
                neighbors, G, k
            )

        return tau

    @staticmethod
    def _compute_k_core_diversity(neighbors, G, k):

        neighbor_list = list(neighbors)
        n_neighbors = len(neighbor_list)
        vertex_to_idx = {neighbor_list[i]: i for i in range(n_neighbors)}

        #adjacency list for neighbor-induced subgraph
        subgraph_adj = [[] for _ in range(n_neighbors)]

        for i, u in enumerate(neighbor_list):
            for w in G.adj_list[u]:
                if w in neighbors:
                    j = vertex_to_idx[w]
                    subgraph_adj[i].append(j)

        #k-core decomposition
        remaining_vertices = kCoreBaseStructuralDiversity._k_core_decomposition(
            subgraph_adj, k
        )

        #connected components in k-core
        return kCoreBaseStructuralDiversity._count_connected_components(
            subgraph_adj, remaining_vertices
        )

    @staticmethod
    def _k_core_decomposition(adj_list, k):
        n = len(adj_list)
        degrees = [len(adj_list[i]) for i in range(n)]
        removed = [False] * n

        #degree < k
        queue = deque()
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)

        #remove vertices with insufficient degree
        while queue:
            i = queue.popleft()
            if removed[i]:
                continue

            removed[i] = True

            #update degrees of neighbors
            for j in adj_list[i]:
                if not removed[j]:
                    degrees[j] -= 1
                    if degrees[j] < k:
                        queue.append(j)

        #return remaining vertices
        return {i for i in range(n) if not removed[i]}

    @staticmethod
    def _count_connected_components(adj_list, vertices):

        if not vertices:
            return 0

        visited = set()
        components = 0

        for start in vertices:
            if start not in visited:
                components += 1

                #DFS
                stack = [start]
                while stack:
                    curr = stack.pop()
                    if curr in visited:
                        continue

                    visited.add(curr)

                    #add unvisited neighbors to stack
                    for neighbor in adj_list[curr]:
                        if neighbor in vertices and neighbor not in visited:
                            stack.append(neighbor)

        return components

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
