#!/usr/bin/env python3
# Auto-generated for 5485390

STUDENT_ID = "5485390"
STUDENT_NAME = "Xiyu Fu"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        for v in range(n):
            nbrs = G.adj_list[v]
            # If the number of neighbors is less than k, there is no k-core
            if len(nbrs) < k:
                sd[v] = 0
                continue
            H_adj = kCoreBaseStructuralDiversity._build_subgraph(G.adj_list, nbrs)
            remaining = kCoreBaseStructuralDiversity._peel_k_core(H_adj, k)
            sd[v] = kCoreBaseStructuralDiversity._count_cc(H_adj, remaining)

        return sd

    @staticmethod
    def _build_subgraph(adj_list, nbrs):
        #Constructs the adjacency structure of the induced subgraph, preserving only the vertices in nbrs and the edges between them.

        nbr_set = set(nbrs)
        H_adj = {u: [] for u in nbrs}
        for u in nbrs:
            for w in adj_list[u]:
                if w in nbr_set:
                    H_adj[u].append(w)
        return H_adj

    @staticmethod
    def _peel_k_core(H_adj, k):
        #Performs k-core stripping on the subgraph H_adj, returning a list of remaining vertices.

        deg = {u: len(neigh) for u, neigh in H_adj.items()}
        queue = deque([u for u, d in deg.items() if d < k])
        removed = set()
        while queue:
            u = queue.popleft()
            if u in removed:
                continue
            removed.add(u)
            for w in H_adj[u]:
                if w not in removed:
                    deg[w] -= 1
                    if deg[w] < k:
                        queue.append(w)
        # Returns undeleted vertices
        return [u for u in H_adj if u not in removed]

    @staticmethod
    def _count_cc(H_adj, remaining):
        #Count the number of connected components of the remaining vertices in the subgraph H_adj.
        visited = set()
        count = 0
        for u in remaining:
            if u not in visited:
                count += 1
                # DFS
                stack = [u]
                visited.add(u)
                while stack:
                    x = stack.pop()
                    for w in H_adj[x]:
                        if w in remaining and w not in visited:
                            visited.add(w)
                            stack.append(w)
        return count


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
