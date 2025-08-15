#!/usr/bin/env python3
# Auto-generated for 5457479

STUDENT_ID = "5457479"
STUDENT_NAME = "Vincent Huang"

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
        tau = [0] * n
        adjacency = G.adj_list

        degG = [len(neigh) for neigh in adjacency]

        for v, neighs in enumerate(adjacency):

            S = {u for u in neighs if degG[u] >= k}
            if not S:
                continue
            tau[v] = kCoreBaseStructuralDiversity._count_k_core_components(adjacency, S, k)

        return tau

    @staticmethod
    def _count_k_core_components(adjacency, S, k):

        sub_adj = {u: [w for w in adjacency[u] if w in S] for u in S}
        deg = {u: len(sub_adj[u]) for u in S}

        queue = deque([u for u, d in deg.items() if d < k])
        removed = set(queue)
        while queue:
            u = queue.popleft()
            for w in sub_adj[u]:
                if w not in removed:
                    deg[w] -= 1
                    if deg[w] < k:
                        removed.add(w)
                        queue.append(w)

        core = S - removed
        if not core:
            return 0

        visited = set()
        components = 0
        for node in core:
            if node not in visited:
                components += 1
                stack = [node]
                visited.add(node)
                while stack:
                    curr = stack.pop()
                    for nbr in sub_adj[curr]:
                        if nbr in core and nbr not in visited:
                            visited.add(nbr)
                            stack.append(nbr)
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
