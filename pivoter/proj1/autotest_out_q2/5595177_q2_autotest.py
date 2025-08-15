#!/usr/bin/env python3
# Auto-generated for 5595177

STUDENT_ID = "5595177"
STUDENT_NAME = "Xinyi Liu"

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
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                sd[v] = 0
                continue

            # Build subgraph G[N(v)]
            subgraph_nodes = set(neighbors)
            subgraph_adj = {u: [] for u in subgraph_nodes}
            node_list = list(subgraph_nodes)
            node_index = {u: idx for idx, u in enumerate(node_list)}
            m = len(node_list)

            # Build adjacency list for the subgraph
            for u in neighbors:
                for nb in G.adj_list[u]:
                    if nb in subgraph_nodes:
                        subgraph_adj[u].append(nb)

            # Compute k-cores in the subgraph
            degrees = [len(subgraph_adj[u]) for u in node_list]
            queue = deque()
            for i in range(m):
                if degrees[i] < k:
                    queue.append(i)

            remaining = m
            is_removed = [False] * m

            while queue and remaining > 0:
                current = queue.popleft()
                if is_removed[current]:
                    continue
                is_removed[current] = True
                remaining -= 1
                u = node_list[current]
                for nb in subgraph_adj[u]:
                    if nb in node_index:
                        idx = node_index[nb]
                        if not is_removed[idx]:
                            degrees[idx] -= 1
                            if degrees[idx] < k:
                                queue.append(idx)

            # Count the number of k-cores (connected components with all degrees >= k)
            if remaining == 0:
                sd[v] = 0
                continue

            # Find connected components in the remaining subgraph
            visited = [False] * m
            core_count = 0
            for i in range(m):
                if not is_removed[i] and not visited[i]:
                    stack = [i]
                    visited[i] = True
                    core_count += 1
                    while stack:
                        current = stack.pop()
                        u = node_list[current]
                        for nb in subgraph_adj[u]:
                            if nb in node_index:
                                idx = node_index[nb]
                                if not is_removed[idx] and not visited[idx]:
                                    visited[idx] = True
                                    stack.append(idx)
            sd[v] = core_count

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
