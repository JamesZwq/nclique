#!/usr/bin/env python3
# Auto-generated for 5567793

STUDENT_ID = "5567793"
STUDENT_NAME = "Mingyuan Liu"

# ======= 学生代码 =======
from collections import deque, defaultdict

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
        neighbors = G.adj_list
        sd = [0] * n  # τ_k(v) for all v

        for v in range(n):
            N_v = neighbors[v]
            if len(N_v) < k:
                continue

            N_set = set(N_v)
            subgraph = {}
            for u in N_v:
                sub_neighbors = [w for w in neighbors[u] if w in N_set and w != u]
                if sub_neighbors:
                    subgraph[u] = set(sub_neighbors)

            if len(subgraph) < k:
                continue

            core_nodes = kCoreBaseStructuralDiversity.extract_k_core(subgraph, k)
            if not core_nodes:
                continue

            visited = set()
            count = 0
            for node in core_nodes:
                if node not in visited:
                    kCoreBaseStructuralDiversity.dfs(node, subgraph, core_nodes, visited)
                    count += 1

            sd[v] = count

        return sd

    @staticmethod
    def extract_k_core(graph, k):
        """
        Removes nodes with degree < k iteratively to extract k-core
        """
        degrees = {u: len(graph[u]) for u in graph}
        queue = deque([u for u, deg in degrees.items() if deg < k])
        alive = set(graph.keys())

        while queue:
            u = queue.popleft()
            alive.discard(u)
            for v in graph[u]:
                if v in alive:
                    degrees[v] -= 1
                    if degrees[v] == k - 1:
                        queue.append(v)

        return alive

    @staticmethod
    def dfs(start, graph, valid_nodes, visited):
        """
        DFS to count connected components in the k-core
        """
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            for neighbor in graph.get(node, []):
                if neighbor in valid_nodes and neighbor not in visited:
                    stack.append(neighbor)

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
