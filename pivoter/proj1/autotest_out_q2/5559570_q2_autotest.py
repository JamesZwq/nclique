#!/usr/bin/env python3
# Auto-generated for 5559570

STUDENT_ID = "5559570"
STUDENT_NAME = "Kai Song"

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
        # """
        # Parameters
        # ----------
        # G : UndirectedUnweightedGraph
        # k : int
        # Returns
        # -------
        # List[int]  # τ_k(v) for all v
        # """
        n = G.vertex_num
        diversity = [0] * n

        for v in range(n):
            neighbours = G.adj_list[v]
            if not neighbours:
                diversity[v] = 0
                continue

            # Construct the neighbor-induced subgraph (excluding the vertex v itself)）
            subgraph = kCoreBaseStructuralDiversity._build_neighbour_subgraph(G, neighbours)


            diversity[v] = kCoreBaseStructuralDiversity._count_k_core(subgraph, k)

        # print(diversity)
        return diversity

    @staticmethod
    def _build_neighbour_subgraph(G, neighbours):
        neighbour_set = set(neighbours)
        subgraph = {u: [] for u in neighbour_set}

        for u in neighbour_set:
            for v in G.adj_list[u]:
                if v in neighbour_set:
                    subgraph[u].append(v)

        return subgraph

    @staticmethod
    def _count_k_core(sub_G, k):
        # Iteratively remove vertices with degree less than k until all remaining vertices have degree ≥ k.

        if not sub_G:
            return 0

        degrees = {u: len(neighbours) for u, neighbours in sub_G.items()}

        # Perform the peeling operation to remove vertices that do not satisfy the k-core condition.
        deleted = kCoreBaseStructuralDiversity._peel_k_core(sub_G, degrees, k)


        remains = [u for u in sub_G if u not in deleted]
        if not remains:
            return 0

        return kCoreBaseStructuralDiversity._count_connected_components(sub_G, deleted)


    @staticmethod
    def _peel_k_core(sub_G, degrees, k):
        # Remove all nodes with degree less than k and update the degrees of their neighbors.
        deleted = set()
        queue = deque()

        for u in sub_G:
            if degrees[u] < k:
                queue.append(u)
                deleted.add(u)

        while queue:
            u = queue.popleft()
            for v in sub_G[u]:
                if v not in deleted:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        deleted.add(v)
                        queue.append(v)

        return deleted


    @staticmethod
    def _count_connected_components(sub_G, deleted):
        # Count the number of connected components in the remaining graph (i.e., nodes that were not deleted)
        visited = set()
        k_core_count = 0

        for node in sub_G:
            if node in deleted or node in visited:
                continue

            k_core_count += 1
            queue = deque()
            queue.append(node)
            visited.add(node)

            while queue:
                u = queue.popleft()
                for v in sub_G[u]:
                    if v not in deleted and v not in visited:
                        visited.add(v)
                        queue.append(v)

        return k_core_count


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
