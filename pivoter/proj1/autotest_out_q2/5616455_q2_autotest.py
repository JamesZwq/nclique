#!/usr/bin/env python3
# Auto-generated for 5616455

STUDENT_ID = "5616455"
STUDENT_NAME = "Yiwei Wu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        result = [0] * n

        for v in range(n):
            nbr = G.adj_list[v]  # Get neighbour list
            if not nbr:
                continue

            # Build neighbour-induced subgraph
            nbr_set = set(nbr)
            subgraph = defaultdict(list)
            for u in nbr:
                for w in G.adj_list[u]:
                    if w in nbr_set:
                        subgraph[u].append(w)

            # Compute degree of each node in subgraph
            degree = {u: len(subgraph[u]) for u in subgraph}

            # Remove nodes with degree < k
            queue = deque([u for u in subgraph if degree[u] < k])
            while queue:
                u = queue.popleft()
                for w in subgraph[u]:
                    if w in degree:
                        degree[w] -= 1
                        if degree[w] < k and w not in queue:
                            queue.append(w)
                if u in degree:
                    del degree[u]

            # Remaining nodes = k-core. Count connected components
            visited = set()
            count = 0
            for u in degree:
                if u not in visited:
                    count += 1
                    q = deque([u])
                    visited.add(u)
                    while q:
                        curr = q.popleft()
                        for nei in subgraph[curr]:
                            if nei in degree and nei not in visited:
                                visited.add(nei)
                                q.append(nei)

            result[v] = count

        return result


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
