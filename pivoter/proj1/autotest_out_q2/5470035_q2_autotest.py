#!/usr/bin/env python3
# Auto-generated for 5470035

STUDENT_ID = "5470035"
STUDENT_NAME = "(Joshua) Tong Zhang"

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
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        n = G.vertex_num
        # 0) Globally deduplicate each adjacency list to remove any parallel edges
        orig_adj = G.adj_list
        adj = [list(dict.fromkeys(neigh)) for neigh in orig_adj]

        sd = [0] * n

        for v in range(n):
            # 1) Locally deduplicate v's neighbor list again
            nbrs = list(dict.fromkeys(adj[v]))
            m = len(nbrs)
            # early exit if there are fewer than k neighbors
            if m < k:
                sd[v] = 0
                continue

            # 2) Build the neighbor-induced subgraph of v
            # Map each neighbor to an index 0..m-1
            idx = {u: i for i, u in enumerate(nbrs)}
            nbrs_set = set(nbrs)
            # sub_adj[i] holds the list of indices j such that (nbrs[i], nbrs[j]) is an edge
            sub_adj = [[] for _ in range(m)]
            for i, u in enumerate(nbrs):
                # Only iterate over neighbors of u that are also in v's neighbor set
                for w in nbrs_set.intersection(adj[u]):
                    sub_adj[i].append(idx[w])
                # Deduplicate this row to ensure no duplicate edges
                sub_adj[i] = list(dict.fromkeys(sub_adj[i]))

            # 3) Perform k-core peeling on the induced subgraph
            # deg[i] = current degree of node i in subgraph
            deg = [len(nei) for nei in sub_adj]
            alive = [True] * m
            # initialize queue with all nodes whose degree < k
            queue = deque(i for i, d in enumerate(deg) if d < k)
            while queue:
                u = queue.popleft()
                if not alive[u]:
                    continue
                alive[u] = False
                # decrement degree of each alive neighbor, and enqueue if it falls below k
                for w in sub_adj[u]:
                    if alive[w]:
                        deg[w] -= 1
                        if deg[w] < k:
                            queue.append(w)

            # 4) Count connected components among the remaining (alive) nodes
            visited = [False] * m
            cnt = 0
            for i in range(m):
                if alive[i] and not visited[i]:
                    cnt += 1
                    stack = [i]
                    visited[i] = True
                    while stack:
                        x = stack.pop()
                        for y in sub_adj[x]:
                            if alive[y] and not visited[y]:
                                visited[y] = True
                                stack.append(y)

            sd[v] = cnt

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
