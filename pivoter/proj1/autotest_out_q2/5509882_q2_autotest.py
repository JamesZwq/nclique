#!/usr/bin/env python3
# Auto-generated for 5509882

STUDENT_ID = "5509882"
STUDENT_NAME = "Dawei Li"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute the k-core-based structural diversity tau_k(v)
        for each vertex v in an undirected, unweighted graph G.
        """
        n = G.vertex_num
        tau = [0] * n

        for v in range(n):
            # build the neighbor-induced subgraph on N(v)
            nbrs = G.adj_list[v]
            if not nbrs:
                # no neighbors => zero components
                tau[v] = 0
                continue

            nbr_set = set(nbrs)
            # sub_adj[u] = list of neighbors of u within nbr_set
            sub_adj = {u: [] for u in nbr_set}
            # deg_map[u] = degree of u in the induced subgraph
            deg_map = {}

            # construct adjacency and degree for H_v
            for u in nbr_set:
                cnt = 0
                for w in G.adj_list[u]:
                    if w in nbr_set:
                        sub_adj[u].append(w)
                        cnt += 1
                deg_map[u] = cnt

            # perform k-core peeling if k > 0
            removed = set()
            if k > 0:
                queue = deque([u for u, d in deg_map.items() if d < k])
                while queue:
                    u = queue.popleft()
                    if u in removed:
                        continue
                    removed.add(u)
                    # decrement degree of its neighbors
                    for w in sub_adj[u]:
                        if w not in removed:
                            deg_map[w] -= 1
                            if deg_map[w] == k - 1:
                                queue.append(w)

            # remaining vertices after peeling
            if k > 0:
                remaining = [u for u in nbr_set if u not in removed]
            else:
                remaining = list(nbr_set)

            # count connected components in the remaining subgraph
            if not remaining:
                tau[v] = 0
            else:
                visited = set()
                comps = 0
                for u in remaining:
                    if u not in visited:
                        comps += 1
                        # BFS/DFS to mark this component
                        stack = [u]
                        visited.add(u)
                        while stack:
                            x = stack.pop()
                            for w in sub_adj[x]:
                                # only traverse nodes still in the k-core
                                if w not in removed and w not in visited:
                                    visited.add(w)
                                    stack.append(w)
                tau[v] = comps

        return tau

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
