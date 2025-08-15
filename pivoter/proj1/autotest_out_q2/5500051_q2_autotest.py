#!/usr/bin/env python3
# Auto-generated for 5500051

STUDENT_ID = "5500051"
STUDENT_NAME = "Yuchen Mei"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################

from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        result = [0] * n

        # calculate the number of k-cores for each vertex
        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                if k >= 1:
                    result[v] = 0
                    continue

            neighbors = set(neighbors)

            # build subgraph
            sub_G = kCoreBaseStructuralDiversity.build_subgraph(G, neighbors)

            # count number of k-cores in the subgraph
            result[v] = kCoreBaseStructuralDiversity.compute_kcore(sub_G, k)

        return result

    @staticmethod
    def build_subgraph(G, neighbors):
        # build subgraph
        sub_G = {}
        for u in neighbors:
            sub_G[u] = []
        for u in neighbors:
            for neighbor in G.adj_list[u]:
                if neighbor in neighbors:
                    sub_G[u].append(neighbor)
        return sub_G

    @staticmethod
    def compute_kcore(sub_G, k):
        # calculate degrees of all nodes in the subgraph
        degrees = {}
        for u in sub_G:
            degrees[u] = len(sub_G[u])

        deleted = set()
        queue = deque()

        # find all nodes with degree < k to start deleting
        for u in sub_G:
            if degrees[u] < k:
                deleted.add(u)
                queue.append(u)

        # propagate deletion, update neighbors' degrees
        while queue:
            u = queue.popleft()
            for v in sub_G[u]:
                if v in deleted:
                    continue
                degrees[v] -= 1
                if degrees[v] < k:
                    deleted.add(v)
                    queue.append(v)

        # collect all remaining nodes
        remain_vertices = []
        for u in sub_G:
            if u not in deleted:
                remain_vertices.append(u)

        # count connected components in the remaining subgraph
        if not remain_vertices:
            return 0
        else:
            return kCoreBaseStructuralDiversity.count_components(sub_G, remain_vertices)

    @staticmethod
    def count_components(sub_G, remain_vertices):
        # count connected components in the remaining subgraph
        count = 0
        visited = set()
        for u in remain_vertices:
            if u not in visited:
                count += 1
                stack = [u]
                visited.add(u)
                while stack:
                    node = stack.pop()
                    for neighbor in sub_G[node]:
                        if neighbor not in visited and neighbor in remain_vertices:
                            visited.add(neighbor)
                            stack.append(neighbor)
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
