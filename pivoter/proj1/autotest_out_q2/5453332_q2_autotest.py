#!/usr/bin/env python3
# Auto-generated for 5453332

STUDENT_ID = "5453332"
STUDENT_NAME = "Chloe Liu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    @staticmethod
    def process(graph, k):
        """
        Parameters
        ----------
        graph : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """

        n = graph.vertex_num
        # in_h: tracks if vertex is in the current neighbor-induced subgraph
        in_h = [False] * n
        # deg_h: stores current degree within the neighbor-induced subgraph
        deg_h = [0] * n
        # timestamp: prevents BFS revisit during same vertex processing
        timestamp = [-1] * n
        # sd: result list storing structural diversity values
        sd = [0] * n

        # process each vertex v in the graph
        for v in range(n):
            # get neighbors of vertex v
            neighbors = graph.adj_list[v]

            # handle isolated vertices (no neighbors)
            if not neighbors:
                sd[v] = 0 # structural diversity is 0 for isolated vertices
                continue

            # STEP 1: mark all neighbors of v as part of the subgraph h
            for u in neighbors:
                in_h[u] = True

            # STEP 2: calculate degrees for each neighbor within h
            for u in neighbors:
                cnt = 0
                # count neighbors of u that are also in h (i.e., in N(v))
                for w in graph.adj_list[u]:
                    if in_h[w]:
                        cnt += 1
                # store degree in induced subgraph
                deg_h[u] = cnt

            # # STEP 3: prune vertices with degree < k using BFS
            q = deque()
            # initial queue setup: add vertices with deg_h < k
            for u in neighbors:
                if in_h[u] and deg_h[u] < k:
                    in_h[u] = False
                    q.append(u)

            # propagate removal effects through the subgraph
            while q:
                u = q.popleft()
                # update neighbors of removed vertex
                for w in graph.adj_list[u]:
                    # only process if w is still in h
                    if in_h[w]:
                        # decrement degree due to u's removal
                        deg_h[w] -= 1
                        # remove if below threshold
                        if deg_h[w] < k:
                            in_h[w] = False
                            q.append(w)

            # STEP 4: count connected components in pruned subgraph
            comp_count = 0
            for u in neighbors:
                # if vertex is in pruned h and not visited in this iteration
                if in_h[u] and timestamp[u] != v:
                    # found new component
                    comp_count += 1
                    queue_comp = deque([u])
                    # mark as visited for this v
                    timestamp[u] = v
                    # # BFS to explore entire connected component
                    while queue_comp:
                        x = queue_comp.popleft()
                        for w in graph.adj_list[x]:
                            # process neighbors in h not yet visited
                            if in_h[w] and timestamp[w] != v:
                                timestamp[w] = v
                                queue_comp.append(w)

            # store component count as structural diversity for v
            sd[v] = comp_count

            # STEP 5: reset in_h for next vertex processing
            for u in neighbors:
                in_h[u] = False

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
