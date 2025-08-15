#!/usr/bin/env python3
# Auto-generated for 5011525

STUDENT_ID = "5011525"
STUDENT_NAME = "Wing Fung Lau"

# ======= 学生代码 =======
'''
τ_k(v), the k-core-based structural diversity, is computed for each node v in the graph.

For each node v:
1. Its neighbors whose degree is less than k are removed.
   This requires repeatedly removing nodes with degree < k. Since we may need to scan all |V| nodes up to |V| times,
   this pruning process for the neighbor-induced subgraph can take up to O(|V|²) time.

2. After pruning, if no neighbors remain, we set τ_k(v) = 0.

3. Otherwise, we build a pruned subgraph induced by the remaining neighbors and compute its connected components
   using BFS. τ_k(v) = the number of connected components.

This process is done for all nodes v ∈ V.
Therefore, the total time complexity is O(|V|³).
'''

from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    def get_core_adj(neighbors, adj_set, k):
            core_adj = {u: adj_set[u].intersection(neighbors) for u in neighbors}

            removed = True
            while removed:
                removed = False
                for u in list(core_adj.keys()):
                    if len(core_adj[u]) < k:
                        for v in core_adj[u]:
                            core_adj[v].discard(u)
                        del core_adj[u]
                        removed = True

            return core_adj


    def bfs(v, visited, adj_list):
        queue = deque([v])
        while queue:
            curr = queue.popleft()
            if curr in visited:
                continue
            visited.add(curr)
            for v in adj_list[curr]:
                if v not in visited:
                    queue.append(v)

    @staticmethod
    def process(G, k):
        cls = kCoreBaseStructuralDiversity
        n = G.vertex_num
        sd = [0] * n
        adj_set = [set(nbrs) for nbrs in G.adj_list]

        for u in range(n):
            neighbors = adj_set[u]
            if not neighbors:
                continue

            core_adj = cls.get_core_adj(neighbors, adj_set, k)
            visited = set()
            count = 0

            for v in core_adj.keys():
                if v not in visited:
                    cls.bfs(v, visited, core_adj)
                    count += 1

            sd[u] = count

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
