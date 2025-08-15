#!/usr/bin/env python3
# Auto-generated for 5559286

STUDENT_ID = "5559286"
STUDENT_NAME = "Yimu Zhu"

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
        def buildSubgraph(neighbors):
            sub_adj = {}
            for u in neighbors:
                sub_adj[u] = set()
            for u in neighbors:
                for v in G.adj_list[u]:
                    if v in sub_adj and v != u:
                        sub_adj[u].add(v)
            return sub_adj

        def k_coreDecom(sub_adj, k):
            degrees = {}
            for node in sub_adj:
                degrees[node] = len(sub_adj[node])
            removed = set()
            changed = True
            while changed:
                changed = False
                for u in list(sub_adj.keys()):
                    if u not in removed and degrees[u] < k:
                        removed.add(u)
                        changed = True
                        for v in sub_adj[u]:
                            if v not in removed:
                                degrees[v] -= 1
            return set(u for u in sub_adj if u not in removed)

        def countComponents(sub_adj, nodes):
            visited = set()
            count = 0
            for u in nodes:
                if u not in visited:
                    count += 1
                    stack = [u]
                    visited.add(u)
                    while stack:
                        curr = stack.pop()
                        for v in sub_adj[curr]:
                            if v in nodes and v not in visited:
                                visited.add(v)
                                stack.append(v)
            return count
        τ = [0] * G.vertex_num
        for v in range(G.vertex_num):
            neighbors = list(G.adj_list[v])
            if not neighbors:
                τ[v] = 0
                continue
            sub_adj = buildSubgraph(neighbors)
            kcore_nodes = k_coreDecom(sub_adj, k)
            if not kcore_nodes:
                τ[v] = 0
            else:
                τ[v] = countComponents(sub_adj, kcore_nodes)
        return τ

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
