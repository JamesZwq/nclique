#!/usr/bin/env python3
# Auto-generated for 5533790

STUDENT_ID = "5533790"
STUDENT_NAME = "Xueying Sun"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Main method to compute τ_k(v) for all nodes v in the graph G.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The input undirected graph, represented using adjacency lists.
        k : int
            The k-core parameter.

        Returns
        -------
        τ : List[int]
            τ_k(v) for each node v (0-based index), where τ_k(v) is the number of
            surviving connected components in the k-core of the neighbor-induced subgraph of v.
        """
        n = G.vertex_num
        τ = [0] * n  # Initialize result list with zeros

        for v in range(n):
            neighbors = G.adj_list[v]  # Step 1: Get neighbors of node v
            if not neighbors:
                continue  # If no neighbors, τ_k(v) = 0 by definition

            # Step 2: Build neighbor-induced subgraph G[N(v)]
            sub_adj = {u: set() for u in neighbors}
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in sub_adj:
                        sub_adj[u].add(w)
                        sub_adj[w].add(u)  # Ensure undirected edge

            # Step 3: Prune the subgraph to get its k-core
            k_core_adj = kCoreBaseStructuralDiversity.k_core_prune(sub_adj, k)

            # Step 4: Count connected components in the pruned k-core
            visited = set()

            def dfs(start_node, component):
                stack = [start_node]
                while stack:
                    node = stack.pop()
                    if node in visited:
                        continue
                    visited.add(node)
                    component.append(node)
                    stack.extend(nei for nei in k_core_adj[node] if nei not in visited)

            count = 0  # Count of surviving components
            for u in k_core_adj:
                if u not in visited:
                    comp = []
                    dfs(u, comp)
                    if comp:
                        count += 1

            τ[v] = count  # Save result τ_k(v)

        return τ

    @staticmethod
    def k_core_prune(adj, k):
        """
        Performs k-core decomposition on the given adjacency list (subgraph).

        Parameters
        ----------
        adj : dict[int, set[int]]
            Adjacency list of the input subgraph.
        k : int
            The core number to prune to.

        Returns
        -------
        adj : dict[int, set[int]]
            The remaining k-core subgraph's adjacency list.
        """
        adj = {u: set(neis) for u, neis in adj.items()}  # Copy to avoid mutation
        deg = {u: len(neis) for u, neis in adj.items()}

        queue = deque([u for u in adj if deg[u] < k])

        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v in adj:
                    adj[v].discard(u)
                    deg[v] -= 1
                    if deg[v] == k - 1:
                        queue.append(v)
            del adj[u]
            del deg[u]

        return adj

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
