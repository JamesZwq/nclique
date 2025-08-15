#!/usr/bin/env python3
# Auto-generated for 5424367

STUDENT_ID = "5424367"
STUDENT_NAME = "Tony Le"

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
        sd = [0] * n

        for v in range(n):
            neighbours = G.adj_list[v]  # Neighbours of the vertex v
            deg_v = len(neighbours)  # Degree of v

            # Handle edge case for a vertex with no neighbours
            if deg_v == 0:
                sd[v] = 0
                continue

            # Build the neighbour-induced subgraph
            # Time complexity: O(deg(v)^2), for each neighbour we check connections to every other neighbour
            neighbour_to_index = {u: i for i, u in enumerate(neighbours)}  # Map to help build induced subgraph
            induced_adj = [[] for _ in range(deg_v)]  # Adjacency list for the induced subgraph
            for i, u in enumerate(neighbours):
                for w in G.adj_list[u]:
                    # Only add edge if w is also a neighbour of v and not a self-loop
                    if w in neighbours and w != u:
                        induced_adj[i].append(neighbour_to_index[w])

            # K-core pruning on the induced subgraph
            # Time complexity: O(deg(v) + |E_v|), where |E_v| is the number of edges in the induced subgraph
            deg = [len(adj) for adj in induced_adj]  # Degree of each vertex in induced subgraph
            in_kcore = [True] * deg_v  # Track whether vertices are still present in k-core
            queue = deque([i for i in range(deg_v) if deg[i] < k])  # Vertices we need to prune
            while queue:
                u = queue.popleft()
                in_kcore[u] = False
                for w in induced_adj[u]:
                    if in_kcore[w]:
                        deg[w] -= 1
                        # If the neighbour degree drops to less than k, we queue for removal
                        if deg[w] == k - 1:
                            queue.append(w)

            # Count connected components in remaining induced subgraph
            # Time complexity: O(deg(v) + |E_v|) for a BFS on the induced subgraph
            visited = [False] * deg_v  # Track visited for BFS traversal
            components = 0
            for i in range(deg_v):
                if in_kcore[i] and not visited[i]:
                    components += 1
                    # BFS to mark all nodes in this component
                    q = deque([i])
                    visited[i] = True
                    while q:
                        curr = q.popleft()
                        for w in induced_adj[curr]:
                            # Only traverse vertices in the k-core and unvisited nodes
                            if in_kcore[w] and not visited[w]:
                                visited[w] = True
                                q.append(w)
            sd[v] = components

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
