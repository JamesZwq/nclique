#!/usr/bin/env python3
# Auto-generated for 5633246

STUDENT_ID = "5633246"
STUDENT_NAME = "Siqi Zhou"

# ======= 学生代码 =======
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        from collections import deque

        cache = {}

        def get_k_cores_in_subgraph(nodes):
            key = frozenset(nodes)
            if key in cache:
                return cache[key]

            index_map = {u: i for i, u in enumerate(nodes)}
            n = len(nodes)

            sub_adj = [[] for _ in range(n)]
            for i, u in enumerate(nodes):
                for v in G.adj_list[u]:
                    if v in index_map:
                        sub_adj[i].append(index_map[v])

            degree = [len(adj) for adj in sub_adj]
            active = [True] * n

            queue = deque(i for i in range(n) if degree[i] < k)
            while queue:
                u = queue.popleft()
                if not active[u]:
                    continue
                active[u] = False
                for v in sub_adj[u]:
                    if active[v]:
                        degree[v] -= 1
                        if degree[v] == k - 1:
                            queue.append(v)

            visited = [False] * n

            def bfs(start):
                q = deque([start])
                visited[start] = True
                while q:
                    u = q.popleft()
                    for v in sub_adj[u]:
                        if active[v] and not visited[v]:
                            visited[v] = True
                            q.append(v)

            count = 0
            for i in range(n):
                if active[i] and not visited[i]:
                    bfs(i)
                    count += 1

            cache[key] = count
            return count

        result = []
        for v in range(G.vertex_num):
            nbr = G.adj_list[v]
            result.append(get_k_cores_in_subgraph(nbr))

        return result

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
