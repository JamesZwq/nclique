#!/usr/bin/env python3
# Auto-generated for 5475430

STUDENT_ID = "5475430"
STUDENT_NAME = "Yang Wang"

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
            Must have G.vertex_num (n) and G.adj_list: List[List[int]]
        k : int
            Which k-core to compute in each neighbour-induced subgraph.
        Returns
        -------
        List[int]  # τ_k(v) for all v=0..n-1
        """
        n = G.vertex_num
        result = [0] * n

        for v in range(n):
            # 1) build neighbour-induced subgraph H on N(v)
            nbrs = G.adj_list[v]
            H_nodes = set(nbrs)
            if not H_nodes:
                result[v] = 0
                continue

            # adjacency for H: only between nodes in H_nodes
            H_adj = {u: set() for u in H_nodes}
            for u in H_nodes:
                for w in G.adj_list[u]:
                    if w in H_nodes:
                        H_adj[u].add(w)

            # 2) peel off all nodes of degree < k to get the k-core of H
            if k > 0:
                queue = deque(u for u in H_nodes if len(H_adj[u]) < k)
                removed = set()
                while queue:
                    u = queue.popleft()
                    if u in removed:
                        continue
                    removed.add(u)
                    # remove u from its neighbours
                    for w in list(H_adj[u]):
                        H_adj[w].remove(u)
                        # if neighbour's degree just dropped below k, schedule it
                        if w not in removed and len(H_adj[w]) < k:
                            queue.append(w)
                    # clear u’s adjacency
                    H_adj[u].clear()
                core_nodes = [u for u in H_nodes if u not in removed]
            else:
                # 0-core is just the full H
                core_nodes = list(H_nodes)

            # 3) count connected components in that remaining subgraph
            visited = set()
            cnt = 0
            for u in core_nodes:
                if u not in visited:
                    cnt += 1
                    dq = deque([u])
                    visited.add(u)
                    while dq:
                        x = dq.popleft()
                        for w in H_adj[x]:
                            if w not in visited:
                                visited.add(w)
                                dq.append(w)
            result[v] = cnt

        return result



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
