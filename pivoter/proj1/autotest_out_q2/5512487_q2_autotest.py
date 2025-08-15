#!/usr/bin/env python3
# Auto-generated for 5512487

STUDENT_ID = "5512487"
STUDENT_NAME = "Mingze Liu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################


class kCoreBaseStructuralDiversity:
    """
    Parameters
    ----------
    G : UndirectedUnweightedGraph
    k : int
    Returns
    -------
    List[int]  # τ_k(v) for all v
    """

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        adj_list = G.adj_list

        group_tag = 0
        neighbour_tags = [-1] * n
        local_indexes = [-1] * n
        seen_tag = [-1] * n

        sd = [0] * n

        for v in range(n):
            nbrs = adj_list[v]
            if not nbrs:
                continue

            group_tag += 1
            kCoreBaseStructuralDiversity._mark(nbrs, group_tag, neighbour_tags, local_indexes)
            degree = kCoreBaseStructuralDiversity._calculate_degrees(nbrs, adj_list, neighbour_tags, group_tag)
            alive = kCoreBaseStructuralDiversity._peel_k_core(nbrs, adj_list, degree, k, neighbour_tags, local_indexes)
            sd[v] = kCoreBaseStructuralDiversity._count_components(nbrs, adj_list, alive, group_tag, neighbour_tags, local_indexes, seen_tag)
            kCoreBaseStructuralDiversity._unmark(nbrs, local_indexes)
        return sd

    @staticmethod
    def _mark(nbrs, tag, neighbour_tags, local_indexes):
        for i, u in enumerate(nbrs):
            neighbour_tags[u] = tag
            local_indexes[u] = i

    @staticmethod
    def _unmark(nbrs, local_indexes):
        for u in nbrs:
            local_indexes[u] = -1

    @staticmethod
    def _calculate_degrees(nbrs, adj_list, neighbour_tags, tag):
        degree = [0] * len(nbrs)
        for i, u in enumerate(nbrs):
            cnt = 0
            for w in adj_list[u]:
                if neighbour_tags[w] == tag:
                    cnt += 1
            degree[i] = cnt
        return degree

    @staticmethod
    def _peel_k_core(nbrs, adj_list, degree, k, neighbour_tags, local_indexes):
        m = len(nbrs)
        alive = [True] * m
        stack = [i for i, d in enumerate(degree) if d < k]

        while stack:
            i = stack.pop()
            if not alive[i]:
                continue
            alive[i] = False
            u = nbrs[i]
            for w in adj_list[u]:
                if neighbour_tags[w] == neighbour_tags[u]:  # same tag
                    j = local_indexes[w]
                    if alive[j]:
                        degree[j] -= 1
                        if degree[j] == k - 1:
                            stack.append(j)
        return alive

    @staticmethod
    def _count_components(nbrs, adj_list, alive, tag, neighbour_tags, local_indexes, seen_tag):
        num = 0
        for i, keep in enumerate(alive):
            u = nbrs[i]
            if not keep or seen_tag[u] == tag:
                continue

            num += 1
            dq = deque([i])
            seen_tag[u] = tag

            while dq:
                cur_i = dq.pop()
                x = nbrs[cur_i]
                for w in adj_list[x]:
                    if neighbour_tags[w] == tag:
                        j = local_indexes[w]
                        if alive[j] and seen_tag[w] != tag:
                            seen_tag[w] = tag
                            dq.append(j)
        return num

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
