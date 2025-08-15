#!/usr/bin/env python3
# Auto-generated for 5500877

STUDENT_ID = "5500877"
STUDENT_NAME = "Haiyu Wang"

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
        sd = [0] * n  # structural diversity array to store τ_k(v) for each vertex v

        for v in range(n):
            neighbors = G.adj_list[v]          # neighbors of vertex v
            subgraph_nodes = set(neighbors)    # nodes in the induced subgraph: the neighbors only

            # Build the induced subgraph adjacency list of neighbors:
            # Only keep edges between neighbors themselves
            subgraph_adj = {u: [] for u in subgraph_nodes}
            for u in subgraph_nodes:
                for w in G.adj_list[u]:
                    if w in subgraph_nodes:
                        subgraph_adj[u].append(w)

            # Perform k-core peeling on the induced subgraph:
            # Initialize degree counts for each node in subgraph
            degree = {u: len(subgraph_adj[u]) for u in subgraph_nodes}
            # Initialize a queue with nodes having degree less than k (to peel off)
            queue = deque([u for u in subgraph_nodes if degree[u] < k])

            # Iteratively remove nodes with degree < k and update degrees of neighbors
            while queue:
                u = queue.popleft()
                for w in subgraph_adj[u]:
                    if w in degree:
                        degree[w] -= 1
                        # If neighbor's degree drops to k-1, it should also be peeled off
                        if degree[w] == k - 1:
                            queue.append(w)
                # Remove the peeled node from degree dictionary
                degree.pop(u)

            # After peeling, degree dictionary contains nodes in the k-core induced subgraph
            # Count connected components in this remaining k-core subgraph
            visited = set()

            def bfs(start):
                q = deque([start])
                visited.add(start)
                while q:
                    u = q.popleft()
                    for w in subgraph_adj[u]:
                        # Visit only nodes still in k-core and not visited yet
                        if w in degree and w not in visited:
                            visited.add(w)
                            q.append(w)

            # For each remaining node, run BFS to find connected components
            for u in degree:
                if u not in visited:
                    bfs(u)
                    sd[v] += 1  # Increment count of connected components for vertex v

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
