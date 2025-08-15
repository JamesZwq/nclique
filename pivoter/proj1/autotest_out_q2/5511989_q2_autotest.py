#!/usr/bin/env python3
# Auto-generated for 5511989

STUDENT_ID = "5511989"
STUDENT_NAME = "Yuchen Wang"

# ======= 学生代码 =======
from collections import defaultdict, deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Parameters:
            G: UndirectedUnweightedGraph Object, including attributes G.vertex_num and G.adj_list
            k: int，k-core Parameter
        Return:
            tau: List[int]，The value τ_k(v) of each vertex v
        """
        n = G.vertex_num
        adj = G.adj_list
        tau = [0] * n

        for v in range(n):
            nbrs = set(adj[v])
            if not nbrs:
                tau[v] = 0
                continue

            # Construct the neighbor-induced subgraph (excluding vertex v itself)
            sub_G = {}
            for u in nbrs:
                sub_G[u] = [w for w in adj[u] if w in nbrs]

            # Call the compute_k_core algorithm
            tau[v] = kCoreBaseStructuralDiversity._compute_k_core(sub_G, k)

        return tau

    @staticmethod
    def _compute_k_core(sub_G, k):
        """
        Input: sub_G represents a graph in adjacency list form, and it has been restricted to the set of neighbors of a certain vertex v.
        Return:  The number of connected components in the k-core
        """
        if not sub_G:
            return 0

        # Initialize the degree of each point
        degrees = {u: len(neigh) for u, neigh in sub_G.items()}
        deleted = set()
        queue = deque([u for u in sub_G if degrees[u] < k])

        # Delete all points with a degree less than k, and update the degrees of their neighboring points.
        while queue:
            u = queue.popleft()
            if u in deleted:
                continue
            deleted.add(u)
            for v in sub_G[u]:
                if v not in deleted:
                    degrees[v] -= 1
                    if degrees[v] == k - 1:
                        queue.append(v)

        # Points in the k-core
        remains = [u for u in sub_G if u not in deleted]
        if not remains:
            return 0

        # Find the number of connected components in the remaining graph (BFS)
        visited = set()
        count = 0
        for node in remains:
            if node in visited:
                continue
            count += 1
            q = deque([node])
            visited.add(node)
            while q:
                u = q.popleft()
                for v in sub_G[u]:
                    if v not in visited and v not in deleted:
                        visited.add(v)
                        q.append(v)

        return count

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
