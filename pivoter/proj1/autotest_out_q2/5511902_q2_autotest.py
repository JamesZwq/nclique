#!/usr/bin/env python3
# Auto-generated for 5511902

STUDENT_ID = "5511902"
STUDENT_NAME = "Chuansheng Xu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def k_core_decomposition(G):
        n = G.vertex_num
        if n == 0:
            return []

        deg = [len(G.adj_list[i]) for i in range(n)]
        max_deg = max(deg) if deg else 0

        bins = [set() for _ in range(max_deg + 1)]
        current_deg = deg.copy()

        # Add the vertex to the corresponding bucket
        for i in range(n):
            d = current_deg[i]
            if d <= max_deg:
                bins[d].add(i)

        core = [0] * n
        idx = 0

        # Process vertices from smallest to largest in degree
        for _ in range(n):
            while idx <= max_deg and not bins[idx]:
                idx += 1
            if idx > max_deg:
                break

            v = bins[idx].pop()
            core[v] = idx

            # Process all neighbors of v
            for w in G.adj_list[v]:
                if current_deg[w] > idx:
                    old_deg = current_deg[w]
                    if old_deg <= max_deg and w in bins[old_deg]:
                        bins[old_deg].remove(w)

                    current_deg[w] -= 1
                    new_deg = current_deg[w]
                    if new_deg <= max_deg:
                        bins[new_deg].add(w)
        return core

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        coreness = kCoreBaseStructuralDiversity.k_core_decomposition(G)

        adj_set = [set(neighbors) for neighbors in G.adj_list]

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                sd[v] = 0
                continue

            # Construct candidate set S
            S = set()
            for u in neighbors:
                if coreness[u] >= k:
                    S.add(u)

            if not S:
                sd[v] = 0
                continue

            # Initialize the degree of subgraph nodes
            deg_in_subgraph = {}
            for u in S:
                common_neighbors = S & adj_set[u]
                deg_in_subgraph[u] = len(common_neighbors)

            # Decomposing k-core
            deleted = set()
            queue = deque()

            for u in S:
                if deg_in_subgraph[u] < k:
                    queue.append(u)
                    deleted.add(u)

            while queue:
                u = queue.popleft()
                for w in adj_set[u]:
                    if w in S and w not in deleted:
                        deg_in_subgraph[w] -= 1
                        if deg_in_subgraph[w] < k:
                            queue.append(w)
                            deleted.add(w)

            # Calculate the number of connected components of the remaining subgraph
            remains = S - deleted
            if not remains:
                sd[v] = 0
                continue

            # BFS
            visited = set()
            count = 0
            for node in remains:
                if node not in visited:
                    count += 1
                    queue_bfs = deque([node])
                    visited.add(node)
                    while queue_bfs:
                        cur = queue_bfs.popleft()
                        for neighbor in adj_set[cur]:
                            if neighbor in remains and neighbor not in visited:
                                visited.add(neighbor)
                                queue_bfs.append(neighbor)
            sd[v] = count

        return sd

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
