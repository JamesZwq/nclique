#!/usr/bin/env python3
# Auto-generated for 5541528

STUDENT_ID = "5541528"
STUDENT_NAME = "Haomiao Xia"

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
            neighbors = G.adj_list[v]
            neighbor_set = set(neighbors)
            if not neighbors:
                sd[v] = 0
                continue

            subgraph = kCoreBaseStructuralDiversity.get_subgraph(G, neighbor_set)

            sd[v] = kCoreBaseStructuralDiversity.count_k_cores(subgraph, k)

        return sd

    @staticmethod
    def get_subgraph(G, neighbor_set):
        subgraph = {}
        for neighbor in neighbor_set:
            subgraph[neighbor] = []
            for v in G.adj_list[neighbor]:
                if v in neighbor_set and v != neighbor:
                    subgraph[neighbor].append(v)
        return subgraph

    @staticmethod
    def count_k_cores(subgraph, k):
        degree = {}
        remove = set()
        q = deque()

        for node in subgraph:
            degree[node] = len(subgraph[node])
            if degree[node] < k:
                q.append(node)
                remove.add(node)

        while q:
            u = q.popleft()
            for v in subgraph[u]:
                degree[v] -= 1
                if degree[v] < k and v not in remove:
                    q.append(v)
                    remove.add(v)

        remain_nodes = []
        for node in subgraph:
            if node not in remove:
                remain_nodes.append(node)

        if remain_nodes:
            k_cores = kCoreBaseStructuralDiversity.get_k_cores(subgraph, remain_nodes, remove)
            return k_cores

        return 0

    @staticmethod
    def get_k_cores(subgraph, remain_nodes, remove):
        count = 0
        visited = set()
        for node in remain_nodes:
            if node not in visited and node not in remove:
                count += 1
                stack = [node]
                visited.add(node)
                while stack:
                    u = stack.pop()
                    for v in subgraph[u]:
                        if v not in visited and v not in remove:
                                stack.append(v)
                                visited.add(v)
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
