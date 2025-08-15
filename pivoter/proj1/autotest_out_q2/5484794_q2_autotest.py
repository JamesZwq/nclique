#!/usr/bin/env python3
# Auto-generated for 5484794

STUDENT_ID = "5484794"
STUDENT_NAME = "Bowen Yao"

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
        adj = G.adj_list
        degree = [len(adj[v]) for v in range(n)]
        removed = [False] * n

        # Global k-core decomposition on the entire graph
        queue = deque([v for v in range(n) if degree[v] < k])
        while queue:
            u = queue.popleft()
            if removed[u]:
                continue
            removed[u] = True
            for v in adj[u]:
                degree[v] -= 1
                if not removed[v] and degree[v] < k:
                    queue.append(v)

        # Initialize result list for tau_k(v)
        sd = [0] * n

        # For each node v, compute tau_k(v)
        for v in range(n):
            if removed[v]:
                continue
            neighbors = [u for u in adj[v] if not removed[u]]
            if not neighbors:
                continue

            # Build the induced subgraph of v's neighbors
            id_map = {u: i for i, u in enumerate(neighbors)}
            m = len(neighbors)
            sub_adj = [[] for _ in range(m)]
            sub_degree = [0] * m

            for u in neighbors:
                uid = id_map[u]
                for w in adj[u]:
                    if w in id_map:
                        wid = id_map[w]
                        sub_adj[uid].append(wid)
                        sub_degree[uid] += 1

            # Perform local k-core decomposition on the neighbor-induced subgraph
            sub_removed = [False] * m
            q = deque([i for i in range(m) if sub_degree[i] < k])
            while q:
                u = q.popleft()
                if sub_removed[u]:
                    continue
                sub_removed[u] = True
                for w in sub_adj[u]:
                    sub_degree[w] -= 1
                    if not sub_removed[w] and sub_degree[w] < k:
                        q.append(w)

            # Count connected components in the remaining subgraph
            visited = [False] * m
            count = 0
            for i in range(m):
                if sub_removed[i] or visited[i]:
                    continue
                q = deque([i])
                visited[i] = True
                while q:
                    u = q.popleft()
                    for w in sub_adj[u]:
                        if not visited[w] and not sub_removed[w]:
                            visited[w] = True
                            q.append(w)
                count += 1

            sd[v] = count

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
