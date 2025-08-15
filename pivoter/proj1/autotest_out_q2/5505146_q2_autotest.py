#!/usr/bin/env python3
# Auto-generated for 5505146

STUDENT_ID = "5505146"
STUDENT_NAME = "Shixun Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        def get_subgraph(neighbors, full):
            sub = defaultdict(set)
            neighbor_set = set(neighbors)
            for u in neighbors:
                for v in full[u]:
                    if v in neighbor_set:
                        sub[u].add(v)
            return sub

        def get_k_core(adj, k):
            degrees = {}
            for u in adj:
                degrees[u] = len(adj[u])

            q = deque()
            for u in adj:
                if degrees[u] < k:
                    q.append(u)
            while q:
                u = q.popleft()
                for v in adj[u]:
                    adj[v].discard(u)
                    if len(adj[v]) == k - 1:
                        q.append(v)
                adj[u] = set()

            result = {}
            for u, neigh in adj.items():
                if len(neigh) >= k:
                    result[u] = neigh
            return result

        def count_cc(adj):
            visited = set()
            count = 0

            def bfs(start):
                queue = deque([start])
                visited.add(start)
                while queue:
                    u = queue.popleft()
                    for v in adj[u]:
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)

            for node in adj:
                if node not in visited:
                    bfs(node)
                    count += 1
            return count

        n = G.vertex_num
        result = []
        for _ in range(n):
            result.append(0)

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                result[v] = 0
                continue
            induced_adj = get_subgraph(neighbors, G.adj_list)
            k_core_adj = get_k_core(induced_adj, k)
            result[v] = count_cc(k_core_adj)

        return result


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
