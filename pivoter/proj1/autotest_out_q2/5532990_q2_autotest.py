#!/usr/bin/env python3
# Auto-generated for 5532990

STUDENT_ID = "5532990"
STUDENT_NAME = "Jiaxu Zhang"

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
        tau_list = [0] * n

        for v in range(n):
            adj_v = G.adj_list[v]
            if len(adj_v) == 0:
                continue

            neighbors = set(adj_v)
            induced_subgraph = kCoreBaseStructuralDiversity._induce_subgraph(G, neighbors)
            tau_list[v] = kCoreBaseStructuralDiversity._count_connected_kcores(induced_subgraph, k)

        return tau_list

    @staticmethod
    def _induce_subgraph(G, node_set):
        # Build the induced subgraph over a given set of vertices.
        induced = {u: [] for u in node_set}
        for u in node_set:
            for v in G.adj_list[u]:
                if v in node_set:
                    induced[u].append(v)
        return induced

    @staticmethod
    def _count_components(graph, valid_nodes, removed):
        # Count number of connected components in the k-core induced subgraph.
        if not valid_nodes:
            return 0

        visited = set()
        components = 0

        for start in valid_nodes:
            if start in visited:
                continue

            components += 1
            queue = deque([start])
            visited.add(start)

            while queue:
                curr = queue.popleft()
                for neighbor in graph[curr]:
                    if neighbor not in visited and neighbor not in removed:
                        visited.add(neighbor)
                        queue.append(neighbor)

        return components

    @staticmethod
    def _count_connected_kcores(graph, k):
        # Perform k-core pruning and count connected k-core components.
        if not graph:
            return 0

        deg = {u: len(adj) for u, adj in graph.items()}
        deleted = set()
        queue = deque(u for u in graph if deg[u] < k)

        deleted.update(queue)

        while queue:
            u = queue.popleft()
            for v in graph[u]:
                if v not in deleted:
                    deg[v] -= 1
                    if deg[v] < k:
                        queue.append(v)
                        deleted.add(v)

        survivors = [u for u in graph if u not in deleted]
        return kCoreBaseStructuralDiversity._count_components(graph, survivors, deleted)

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
