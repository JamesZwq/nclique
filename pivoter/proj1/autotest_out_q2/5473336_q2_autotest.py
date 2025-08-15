#!/usr/bin/env python3
# Auto-generated for 5473336

STUDENT_ID = "5473336"
STUDENT_NAME = "Junzhe Ren"

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
        n = G.vertex_num

        # Step 1: Compute the global k-core of the graph
        # initial degrees
        degree = [len(G.adj_list[v]) for v in range(n)]
        # track if a node remains in the k-core
        in_kcore = [True] * n

        queue = deque()
        for v in range(n):
            if degree[v] < k:
                queue.append(v)
                # mark for remove
                in_kcore[v] = False

        # Iteratively remove vertuces with degree < k
        while queue:
            v = queue.popleft()
            for u in G.adj_list[v]:
                if in_kcore[u]:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        in_kcore[u] = False

        # Step 2: Compute τ_k(v) for each vertex v
        # final result array
        sd = [0] * n
        for v in range(n):
            if not in_kcore[v]:
                # if v itself is not in the k-core, τ_k(v) = 0
                sd[v] = 0
                continue

            # collect v's neighbors that are also in the k-core
            neighbors = [u for u in G.adj_list[v] if in_kcore[u]]
            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # build induced subgraph over v's neighbors
            neighbor_set = set(neighbors)
            neighbor_to_idx = {u: i for i, u in enumerate(neighbors)}
            subgraph_adj = [[] for _ in range(len(neighbors))]

            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbor_set and w != v:
                        # add edge in subgraph
                        subgraph_adj[neighbor_to_idx[u]].append(neighbor_to_idx[w])

            # Step 2.1: Compute the k-core of the neighbor-induced subgraph
            sub_n = len(neighbors)
            sub_degree = [len(subgraph_adj[i]) for i in range(sub_n)]
            sub_in_kcore = [True] * sub_n

            sub_queue = deque()
            for i in range(sub_n):
                if sub_degree[i] < k:
                    sub_queue.append(i)
                    sub_in_kcore[i] = False

            while sub_queue:
                i = sub_queue.popleft()
                for j in subgraph_adj[i]:
                    if sub_in_kcore[j]:
                        sub_degree[j] -= 1
                        if sub_degree[j] < k:
                            sub_queue.append(j)
                            sub_in_kcore[j] = False

            # Step 2.2: Count number of connected components in the subgraph k-core
            visited = [False] * sub_n
            num_k_cores = 0

            for i in range(sub_n):
                if sub_in_kcore[i] and not visited[i]:
                    # Start a BFS to count this component
                    num_k_cores += 1
                    bfs_queue = deque([i])
                    visited[i] = True
                    while bfs_queue:
                        curr = bfs_queue.popleft()
                        for next_node in subgraph_adj[curr]:
                            if sub_in_kcore[next_node] and not visited[next_node]:
                                visited[next_node] = True
                                bfs_queue.append(next_node)

            # τ_k(v) = number of connected k-core components among v's neighbors
            sd[v] = num_k_cores

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
