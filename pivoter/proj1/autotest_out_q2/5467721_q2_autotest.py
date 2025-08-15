#!/usr/bin/env python3
# Auto-generated for 5467721

STUDENT_ID = "5467721"
STUDENT_NAME = "Zhongwei Kou"

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
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if len(neighbors) == 0:
                sd[v] = 0
                continue

            vertex_map = {u: i for i, u in enumerate(neighbors)}
            sub_n = len(neighbors)

            # Build adjacency list for neighbor-induced subgraph
            sub_adj = [[] for _ in range(sub_n)]
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in vertex_map and w != v:
                        sub_adj[vertex_map[u]].append(vertex_map[w])

            # Find k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity.find_k_cores(sub_adj, k)
            sd[v] = k_cores

        return sd

    @staticmethod
    def find_k_cores(adj_list, k):
        """
        Find the number of k-cores in a graph
        """
        n = len(adj_list)
        if n == 0:
            return 0

        # Special case for k=0
        if k == 0:
            visited = [False] * n
            component_count = 0

            for i in range(n):
                if not visited[i]:
                    queue = deque([i])
                    visited[i] = True

                    while queue:
                        u = queue.popleft()
                        for v in adj_list[u]:
                            if not visited[v]:
                                visited[v] = True
                                queue.append(v)

                    component_count += 1

            return component_count

        degree = [len(adj_list[i]) for i in range(n)]
        active = [True] * n

        # Iteratively remove vertices with degree < k
        queue = deque()
        for i in range(n):
            if degree[i] < k:
                queue.append(i)
                active[i] = False  # Mark as inactive immediately

        while queue:
            u = queue.popleft()

            for v in adj_list[u]:
                if active[v]:
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)
                        active[v] = False

        #Ensure all active vertices have degree >= k in the remaining subgraph
        for i in range(n):
            if active[i]:
                actual_degree = sum(1 for j in adj_list[i] if active[j])
                if actual_degree < k:
                    active[i] = False

        # Count connected components in remaining graph (k-core)
        visited = [False] * n
        k_core_count = 0

        for i in range(n):
            if active[i] and not visited[i]:
                kCoreBaseStructuralDiversity.bfs(i, adj_list, active, visited)
                k_core_count += 1

        return k_core_count

    @staticmethod
    def bfs(start, adj_list, active, visited):
        """
        BFS to mark all vertices in a connected component
        """
        queue = deque([start])
        visited[start] = True

        while queue:
            u = queue.popleft()
            for v in adj_list[u]:
                if active[v] and not visited[v]:
                    visited[v] = True
                    queue.append(v)

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
