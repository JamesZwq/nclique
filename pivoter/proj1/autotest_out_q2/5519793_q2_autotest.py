#!/usr/bin/env python3
# Auto-generated for 5519793

STUDENT_ID = "5519793"
STUDENT_NAME = "Yingxin Li"

# ======= 学生代码 =======
from collections import deque


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
            nbrs = G.adj_list[v]  # adj_list下标即顶点序号
            if not nbrs:
                sd[v] = 0
                continue

            # 只保留 v 的邻居之间的边
            nbr_set = set(nbrs)   # 变成集合可以优化性能，加速
            # 找该点邻居集合范围中的每个点的邻居子集
            sub_adj = {u: [] for u in nbr_set}
            for u in nbr_set:
                for w in G.adj_list[u]:
                    if w in nbr_set:
                        sub_adj[u].append(w)

            # 记录每个点在子集中的度数
            deg = {u: len(sub_adj[u]) for u in nbr_set}
            # 先把小于k的提取出来，然后先删除再继续看，不能直接提取大于k的
            q = deque()
            alive = set(nbr_set)
            for u in nbr_set:
                deg[u] = len(sub_adj[u])
                if deg[u] < k:
                    q.append(u)

            while q:
                u = q.popleft()
                if u not in alive:
                    continue
                alive.discard(u)
                for w in sub_adj[u]:
                    if w not in alive:
                        continue
                    deg[w] -= 1
                    if deg[w] < k:
                        q.append(w)

            # 剩余节点
            visited = set()
            comp_num = 0
            q_comp = deque()

            for v0 in alive:
                if v0 in visited:
                    continue
                visited.add(v0)
                comp_num += 1
                q_comp.append(v0)
                # bfs
                while q_comp:
                    x = q_comp.popleft()
                    for y in sub_adj[x]:
                        if y in alive and y not in visited:
                            visited.add(y)
                            q_comp.append(y)

            sd[v] = comp_num


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
