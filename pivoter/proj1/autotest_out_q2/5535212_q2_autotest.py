#!/usr/bin/env python3
# Auto-generated for 5535212

STUDENT_ID = "5535212"
STUDENT_NAME = "Mingrunlin Wang"

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
        # For each v, build its neighbor subgraph, peel its k-core, then count components
        for v in range(n):
            N = set(G.adj_list[v])
            if not N or (k > 0 and len(N) <= k):
                # no neighbors or too few to form any k-core
                sd[v] = 0
                continue

            # compute degree of each u ∈ N within the subgraph induced by N
            deg = {u: 0 for u in N}
            for u in N:
                for w in G.adj_list[u]:
                    if w in deg:
                        deg[u] += 1

            # peel off vertices with deg < k
            queue = deque(u for u, d in deg.items() if d < k)
            removed = set()
            while queue:
                u = queue.popleft()
                if u in removed:
                    continue
                removed.add(u)
                for w in G.adj_list[u]:
                    if w in deg and w not in removed:
                        deg[w] -= 1
                        if deg[w] == k - 1:
                            queue.append(w)

            # survivors = the k-core of the neighbor subgraph
            survivors = [u for u in N if u not in removed]
            if not survivors:
                sd[v] = 0
                continue

            # count connected components among survivors
            visited = set()
            comp_count = 0
            for u in survivors:
                if u not in visited:
                    comp_count += 1
                    bfs = deque([u])
                    visited.add(u)
                    while bfs:
                        cur = bfs.popleft()
                        for w in G.adj_list[cur]:
                            if w in N and w not in removed and w not in visited:
                                visited.add(w)
                                bfs.append(w)

            sd[v] = comp_count

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
