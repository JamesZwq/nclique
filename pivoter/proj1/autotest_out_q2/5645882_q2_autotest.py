#!/usr/bin/env python3
# Auto-generated for 5645882

STUDENT_ID = "5645882"
STUDENT_NAME = "Jianing Liu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass
    @staticmethod
    def process(G, k):
        # Prepare a list of neighbor‐sets
        neighbor_sets = [set(adj) for adj in G.adj_list]
        n = G.vertex_num

        # hold the number of k‐cores in v's neighbor‐induced subgraph
        tau = [0] * n

        for v in range(n):
            nbrs = neighbor_sets[v]
            if not nbrs:
                tau[v] = 0
                continue

            # build H, the subgraph induced by v's neighbors
            H = {}
            for u in nbrs:
                H[u] = neighbor_sets[u].intersection(nbrs)

            # perform k‐core peeling on H
            degree = {u: len(adj) for u, adj in H.items()}
            queue = deque(u for u, d in degree.items() if d < k)
            removed = set(queue)

            while queue:
                cur = queue.popleft()
                for nei in H[cur]:
                    if nei not in removed:
                        degree[nei] -= 1
                        if degree[nei] < k:
                            removed.add(nei)
                            queue.append(nei)

            core_nodes = set(H.keys()) - removed

            # count connected components among the core nodes
            seen = set()
            count_cc = 0
            for u in core_nodes:
                if u not in seen:
                    count_cc += 1
                    dq = deque([u])
                    seen.add(u)
                    while dq:
                        x = dq.popleft()
                        for y in H[x]:
                            if y in core_nodes and y not in seen:
                                seen.add(y)
                                dq.append(y)

            tau[v] = count_cc

        return tau

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
