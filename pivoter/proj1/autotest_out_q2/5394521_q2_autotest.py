#!/usr/bin/env python3
# Auto-generated for 5394521

STUDENT_ID = "5394521"
STUDENT_NAME = "Tianbo Liu"

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
        adj = G.adj_list
        tau = [0] * n

        for v in range(n):

            Nv_set = set(adj[v])
            if v in Nv_set:
                Nv_set.remove(v)
            m = len(Nv_set)
            if k > 0 and m < k + 1:
                continue

            Nv = list(Nv_set)
            idx = {u: i for i, u in enumerate(Nv)}

            g_set = [set() for _ in range(m)]
            for i, u in enumerate(Nv):
                for w in adj[u]:
                    j = idx.get(w)
                    if j is not None and j != i:
                        g_set[i].add(j)
                        g_set[j].add(i)

            g   = [list(s) for s in g_set]
            deg = [len(s) for s in g_set]

            removed = [False] * m
            if k > 0:
                q = deque(range(m))
                while q:
                    i = q.popleft()
                    if removed[i] or deg[i] >= k:
                        continue
                    removed[i] = True
                    for j in g[i]:
                        deg[j] -= 1
                        if not removed[j]:
                            q.append(j)

            comp_cnt = 0
            visited  = [False] * m
            for i in range(m):
                if removed[i] or visited[i]:
                    continue
                comp_cnt += 1
                stack = [i]; visited[i] = True
                while stack:
                    s = stack.pop()
                    for j in g[s]:
                        if not removed[j] and not visited[j]:
                            visited[j] = True
                            stack.append(j)

            tau[v] = comp_cnt

        return tau


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
