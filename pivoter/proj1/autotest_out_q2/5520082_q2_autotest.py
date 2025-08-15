#!/usr/bin/env python3
# Auto-generated for 5520082

STUDENT_ID = "5520082"
STUDENT_NAME = "Yian Zhu"

# ======= 学生代码 =======
from collections import deque

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

            neighbor_set = set(neighbors)
            subgraph_nodes = neighbors

            degrees = {}
            node_to_index = {u: i for i, u in enumerate(subgraph_nodes)}
            index_to_node = {i: u for i, u in enumerate(subgraph_nodes)}
            subgraph_size = len(subgraph_nodes)

            for u in subgraph_nodes:
                degrees[u] = 0

            for u in subgraph_nodes:
                for neighbor in G.adj_list[u]:
                    if neighbor in neighbor_set:
                        degrees[u] += 1
            q = deque()
            for u in subgraph_nodes:
                if degrees[u] < k:
                    q.append(u)
            while q:
                u = q.popleft()
                if u not in degrees:
                    continue
                del degrees[u]
                for neighbor in G.adj_list[u]:
                    if neighbor in degrees:
                        degrees[neighbor] -= 1
                        if degrees[neighbor] < k and neighbor not in q:
                            q.append(neighbor)

            remaining_nodes = set(degrees.keys())
            visited = set()
            components = 0

            for u in subgraph_nodes:
                if u in remaining_nodes and u not in visited:
                    components += 1
                    queue = deque([u])
                    visited.add(u)

                    while queue:
                        current = queue.popleft()
                        for neighbor in G.adj_list[current]:
                            if neighbor in remaining_nodes and neighbor not in visited:
                                visited.add(neighbor)
                                queue.append(neighbor)

            sd[v] = components

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
