#!/usr/bin/env python3
# Auto-generated for 5536400

STUDENT_ID = "5536400"
STUDENT_NAME = "Yufei Zhang"

# ======= 学生代码 =======
from collections import deque, defaultdict

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Computes τ_k(v) for every vertex v in the graph G.
        τ_k(v) is the number of connected components in the k-core
        of the subgraph induced by v's neighbors.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The input undirected, unweighted graph with an adjacency list.
        k : int
            The core value k for k-core-based structural diversity.

        Returns
        -------
        List[int]
            A list τ where τ[v] = τ_k(v) for all nodes v in G.
        """

        n = G.vertex_num  # Total number of nodes
        result = [0] * n  # Final τ_k(v) for each node v

        # Process each node v independently
        for v in range(n):

            # Step 1: Get the neighbors of v
            nbrs = G.adj_list[v]
            if not nbrs:
                result[v] = 0
                continue

            # Step 2: Build the neighbor-induced subgraph of v
            # Only include edges between v's neighbors
            induced_adj = defaultdict(set)
            nbr_set = set(nbrs)
            for u in nbrs:
                for w in G.adj_list[u]:
                    if w in nbr_set:
                        induced_adj[u].add(w)
                        induced_adj[w].add(u)

            # Step 3: Perform k-core peeling on the induced subgraph
            # Initialize degrees and active node status
            degrees = {u: len(induced_adj[u]) for u in nbr_set}
            in_kcore = {u: True for u in nbr_set}
            queue = deque()

            # Initial removal of nodes with degree < k
            for u in nbr_set:
                if degrees[u] < k:
                    queue.append(u)
                    in_kcore[u] = False

            # Iteratively peel nodes with degree < k
            while queue:
                u = queue.popleft()
                for w in induced_adj[u]:
                    if in_kcore.get(w, False):
                        degrees[w] -= 1
                        if degrees[w] < k:
                            in_kcore[w] = False
                            queue.append(w)

            # Step 4: Count connected components in the k-core subgraph
            visited = set()
            count = 0  # Number of connected components

            for u in nbr_set:
                if in_kcore.get(u, False) and u not in visited:
                    count += 1
                    # BFS to traverse the current component
                    q = deque()
                    q.append(u)
                    visited.add(u)

                    while q:
                        curr = q.popleft()
                        for nb in induced_adj[curr]:
                            if in_kcore.get(nb, False) and nb not in visited:
                                visited.add(nb)
                                q.append(nb)

            # Step 5: Set τ_k(v) = number of connected components in k-core
            result[v] = count

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
