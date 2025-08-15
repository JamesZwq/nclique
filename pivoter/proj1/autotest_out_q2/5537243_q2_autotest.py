#!/usr/bin/env python3
# Auto-generated for 5537243

STUDENT_ID = "5537243"
STUDENT_NAME = "Dongshuai Ding"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    @staticmethod
    def process(G, k):

        #Compute, for each vertex v in graph G, the number of k-core connected components in its neighbor-induced subgraph.
        n = G.vertex_num
        adj = G.adj_list
        cache = {}
        result = [0] * n

        # Iterate over all vertices and compute the number of k-core connected components
        for v in range(n):
            nbrs = frozenset(adj[v])
            # Prune: if a vertex has fewer than k neighbors, it cannot form a k-core
            if len(nbrs) < k:
                result[v] = 0
                continue
            if nbrs in cache:
                result[v] = cache[nbrs]
                continue

            # Construct the neighbor-induced subgraph
            subG = kCoreBaseStructuralDiversity._build_induced_graph(adj, nbrs)
            # Perform k-core decomposition on the subgraph and count the connected components
            count = kCoreBaseStructuralDiversity._k_core_components(subG, k)
            cache[nbrs] = count
            result[v] = count

        return result

    @staticmethod
    def _build_induced_graph(adj_list, nodes):
        # Construct the neighbor-induced subgraph: For each u in nodes, filter its adjacency list to include only neighbors that are also in nodes
        return {
            u: [w for w in adj_list[u] if w in nodes]
            for u in nodes
        }

    @staticmethod
    def _k_core_components(subG, k):

        # Perform k-core peeling on subG and count the number of connected components formed by the remaining nodes
        # Initialize the degree of each node
        degrees = {u: len(neighbors) for u, neighbors in subG.items()}
        # Enqueue nodes with degree < k and mark them as removed
        queue = deque([u for u, deg in degrees.items() if deg < k])
        removed = set(queue)

        # Iteratively peel: for each removed node u, update the degrees of its neighbors
        while queue:
            u = queue.popleft()
            for w in subG[u]:
                if w not in removed:
                    degrees[w] -= 1
                    if degrees[w] < k:
                        removed.add(w)
                        queue.append(w)

        # The remaining nodes constitute the k-core
        core_nodes = [u for u in subG if u not in removed]
        if not core_nodes:
            return 0

        # Count the number of connected components among the remaining nodes using BFS
        comp_count = 0
        visited = set()
        for u in core_nodes:
            if u in visited:
                continue
            comp_count += 1
            bfs = deque([u])
            visited.add(u)
            while bfs:
                x = bfs.popleft()
                for w in subG[x]:
                    if w not in removed and w not in visited:
                        visited.add(w)
                        bfs.append(w)
        return comp_count

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
