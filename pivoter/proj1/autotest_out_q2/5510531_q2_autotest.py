#!/usr/bin/env python3
# Auto-generated for 5510531

STUDENT_ID = "5510531"
STUDENT_NAME = "Yuqi Qiu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from collections import defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0 for _ in range(n)]

        for v in range(n):
            neighbors = G.adj_list[v]
            neighbor_count = len(neighbors) if neighbors else 0
            if neighbor_count == 0 or neighbor_count < k:
                sd[v] = 0
                continue

            neighbor_set = set(neighbors)
            subgraph = defaultdict(list)
            for u in neighbor_set:
                if u >= n:
                    continue
                current_neighbors = G.adj_list[u]
                if not current_neighbors:
                    continue
                for w in current_neighbors:
                    if w in neighbor_set and u < w:
                        subgraph[u].append(w)
                        subgraph[w].append(u)

            if len(subgraph) == 0:
                sd[v] = 0
                continue

            degree = {}
            for node in subgraph:
                degree[node] = len(subgraph[node])
            to_delete = deque([node for node in subgraph if degree[node] < k])

            while len(to_delete) > 0:
                current_node = to_delete.popleft()
                for neighbor in subgraph.get(current_node, []):
                    if neighbor in subgraph:
                        subgraph[neighbor].remove(current_node)
                        degree[neighbor] -= 1
                        if degree[neighbor] == k - 1:
                            to_delete.append(neighbor)
                if current_node in subgraph:
                    del subgraph[current_node]
                if current_node in degree:
                    del degree[current_node]

            if len(subgraph) == 0:
                sd[v] = 0
                continue

            visited = set()
            components = 0
            for node in subgraph:
                if node in visited:
                    continue
                components += 1
                queue = deque()
                queue.append(node)
                visited.add(node)

                while len(queue) > 0:
                    current = queue.popleft()
                    neighbors_current = subgraph[current]
                    if not neighbors_current:
                        continue
                    for neighbor in neighbors_current:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

            sd[v] = components

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
