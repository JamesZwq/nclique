#!/usr/bin/env python3
# Auto-generated for 5474114

STUDENT_ID = "5474114"
STUDENT_NAME = "Letong Hu"

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

        tau_values = [0] * G.vertex_num

        for v in range(G.vertex_num):
            nbrs = G.adj_list[v]

            # If no neighbors, diversity is zero
            if not nbrs:
                continue

            # Build induced subgraph using neighbors of v
            induced_graph = {node: set() for node in nbrs}
            for node in nbrs:
                induced_graph[node].update(w for w in G.adj_list[node] if w in induced_graph)

            # Perform k-core pruning
            # Use a queue for nodes that might need removal
            prune_queue = deque([n for n in induced_graph if len(induced_graph[n]) < k])
            while prune_queue:
                n = prune_queue.popleft()
                if n not in induced_graph:
                    continue
                for adj in induced_graph[n]:
                    induced_graph[adj].discard(n)
                    if len(induced_graph[adj]) == k-1:
                        prune_queue.append(adj)
                del induced_graph[n]

            # Count connected components in the remaining graph
            seen, comp_count = set(), 0
            for node in induced_graph:
                if node in seen:
                    continue
                comp_count += 1
                q = deque([node])
                while q:
                    cur = q.pop()
                    if cur in seen:
                        continue
                    seen.add(cur)
                    q.extend(induced_graph[cur] - seen)

            tau_values[v] = comp_count

        return tau_values

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
