#!/usr/bin/env python3
# Auto-generated for 5471987

STUDENT_ID = "5471987"
STUDENT_NAME = "Bohua Wang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from typing import List, Dict   # add
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
        
        n = G.vertex_num          # node count
        tau = [0] * n             # result

        v = 0                     
        while v < n:
            neigh_list = G.adj_list[v]    # list of neighbours

            # no neighbour → tau[v] stays 0
            if len(neigh_list) == 0:
                v += 1
                continue

            # put neighbours into a set
            neigh_set = set()
            i = 0
            while i < len(neigh_list):
                neigh_set.add(neigh_list[i])
                i += 1

            # make empty sub-graph dict
            sub_adj = {}
            for node in neigh_set:
                sub_adj[node] = []        # empty list first

            # fill sub-graph edges
            for node in neigh_set:
                full_edges = G.adj_list[node]
                j = 0
                while j < len(full_edges):
                    other = full_edges[j]
                    if other in neigh_set:
                        sub_adj[node].append(other)
                    j += 1

            # count parts inside k-core
            parts = kCoreBaseStructuralDiversity._count_parts(sub_adj, k)
            tau[v] = parts

            v += 1                        # next v

        return tau

    # helper: count components in k-core
    @staticmethod
    def _count_parts(sub_adj, k):
        # empty sub-graph
        if len(sub_adj) == 0:
            return 0

        # degree table
        degree = {}
        for node in sub_adj:
            degree[node] = len(sub_adj[node])

        # list of nodes to delete
        to_del = []
        for node in degree:
            if degree[node] < k:
                to_del.append(node)

        removed = set()

        # peel nodes
        while len(to_del) > 0:
            cur = to_del.pop()            # last one
            if cur in removed:
                continue
            removed.add(cur)

            # touch every neighbour
            idx = 0
            while idx < len(sub_adj[cur]):
                nb = sub_adj[cur][idx]
                idx += 1
                if nb in removed:
                    continue
                degree[nb] -= 1
                if degree[nb] < k:
                    to_del.append(nb)

        # collect nodes still alive
        left = []
        for node in sub_adj:
            if node not in removed:
                left.append(node)

        if len(left) == 0:
            return 0

        # BFS to count parts
        seen = set()
        parts = 0
        left_idx = 0
        while left_idx < len(left):
            start = left[left_idx]
            left_idx += 1
            if start in seen:
                continue

            parts += 1                    # new part
            q = deque()
            q.append(start)
            seen.add(start)

            while len(q) > 0:
                head = q.popleft()
                for nb in sub_adj[head]:
                    if nb in removed:
                        continue
                    if nb in seen:
                        continue
                    seen.add(nb)
                    q.append(nb)

        return parts
    
    
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
