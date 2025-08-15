#!/usr/bin/env python3
# Auto-generated for 5502828

STUDENT_ID = "5502828"
STUDENT_NAME = "Shuobo Wang"

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
            # Get the neighbor set of vertex v
            neighbors = set()
            adj_v = G.adj_list[v]
            for u in adj_v:
                neighbors.add(u)

            # Skip if v has no neighbors
            if len(neighbors) == 0:
                continue

            # Build the adjacency list of the neighbor-induced subgraph
            subgraph = {}
            for u in neighbors:
                subgraph[u] = []

            for u in neighbors:
                adj_u = G.adj_list[u]
                for w in adj_u:
                    if w in neighbors:
                        subgraph[u].append(w)

            # Compute k-core-based structural diversity
            sd[v] = kCoreBaseStructuralDiversity.count_k_core(subgraph, k)

        return sd

    @staticmethod
    def count_k_core(graph, k):
        if not graph:
            return 0

        # Initialize degree dictionary
        deg = {}
        for u in graph:
            neighbors = graph[u]
            deg[u] = 0
            for _ in neighbors:
                deg[u] += 1

        # Find all nodes with degree < k
        removed = set()
        queue = deque()
        for u in deg:
            if deg[u] < k:
                removed.add(u)
                queue.append(u)

        # Remove nodes iteratively
        while len(queue) > 0:
            u = queue.popleft()
            for v in graph[u]:
                if v not in removed:
                    deg[v] -= 1
                    if deg[v] < k:
                        removed.add(v)
                        queue.append(v)

        # Filter remaining nodes
        remaining = []
        for u in graph:
            if u not in removed:
                remaining.append(u)

        if len(remaining) == 0:
            return 0

        # Count the number of connected components in remaining nodes
        visited = set()
        count = 0
        for u in remaining:
            if u in visited:
                continue
            count += 1
            queue = deque()
            queue.append(u)
            visited.add(u)

            while len(queue) > 0:
                curr = queue.popleft()
                for nei in graph[curr]:
                    if nei not in removed and nei not in visited:
                        visited.add(nei)
                        queue.append(nei)

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
