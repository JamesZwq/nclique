#!/usr/bin/env python3
# Auto-generated for 5459790

STUDENT_ID = "5459790"
STUDENT_NAME = "Hao Sun"

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
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                sd[v] = 0
                continue

            # 构造邻居诱导子图
            subgraph = {}
            for u in neighbors:
                subgraph[u] = [w for w in G.adj_list[u] if w in neighbors]

            # 计算 k-core
            k_core_graph = kCoreBaseStructuralDiversity.compute_k_core(subgraph, k)

            # 统计 k-core 中连通分量个数
            num_components = kCoreBaseStructuralDiversity.count_connected_components(k_core_graph)

            sd[v] = num_components

        return sd

    @staticmethod
    @staticmethod
    def compute_k_core(adj_list, k):
        adj_copy = {u: list(neigh) for u, neigh in adj_list.items()}
        degrees = {u: len(adj_copy[u]) for u in adj_copy}
        queue = deque([u for u in adj_copy if degrees[u] < k])
        while queue:
            u = queue.popleft()
            for v in adj_copy[u]:
                if v in adj_copy:
                    degrees[v] -= 1
                    if degrees[v] == k-1:
                        queue.append(v)
            adj_copy.pop(u)
        return adj_copy


    @staticmethod
    @staticmethod
    @staticmethod
    def count_connected_components(adj_list):
        visited = set()
        count = 0
        for u in adj_list:
            if u not in visited:
                queue = deque([u])
                visited.add(u)
                while queue:
                    v = queue.popleft()
                    for w in adj_list.get(v, []):
                        if w in adj_list and w not in visited:
                            visited.add(w)
                            queue.append(w)
                count += 1
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
