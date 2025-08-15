#!/usr/bin/env python3
# Auto-generated for 5489901

STUDENT_ID = "5489901"
STUDENT_NAME = "Ziyi Shi"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    def _process_single_vertex(args):
        # Unpack the arguments: adjacency list (G_adj), vertex index (v), and k value
        G_adj, v, k = args
        neighbors = G_adj[v]  # List of neighbors for vertex v
        d = len(neighbors)    # Degree of vertex v (number of neighbors)
        # If the degree is less than k or has no neighbors, this vertex cannot have a k-core among its neighbors
        if d < k or d == 0:
            return 0

        # Build a mapping from each neighbor node to its index in the induced subgraph
        sub_id_map = {u: i for i, u in enumerate(neighbors)}
        # Create adjacency list for the induced subgraph of neighbors
        induced_adj = [set() for _ in range(d)]
        # Store the degree of each node in the induced subgraph
        deg = [0] * d
        for i, u in enumerate(neighbors):
            for w in G_adj[u]:
                # Only consider neighbors that are also in the subgraph and not itself
                if w in sub_id_map and w != u:
                    induced_adj[i].add(sub_id_map[w])
            deg[i] = len(induced_adj[i])  # Degree in the induced subgraph

        # Pruning: If all neighbors have degree less than k in the induced subgraph, no k-core exists
        if all(deg[i] < k for i in range(d)):
            return 0

        # alive[i] indicates whether node i in the induced subgraph is still in the k-core
        alive = [True] * d
        # Initialize queue with nodes whose degree is less than k (to be removed)
        q = deque(i for i in range(d) if deg[i] < k)
        for i in q:
            alive[i] = False  # Mark as removed
        removed = len(q)     # Count of removed nodes

        # Perform k-core peeling: iteratively remove nodes with degree < k
        while q:
            u = q.popleft()
            for w in induced_adj[u]:
                if alive[w]:
                    deg[w] -= 1  # Decrease degree of neighbor
                    if deg[w] < k:
                        alive[w] = False
                        q.append(w)
                        removed += 1
        # Pruning: If all nodes are removed, there is no k-core
        if removed == d:
            return 0

        # Count the number of connected components among the remaining nodes in the k-core
        visited = [False] * d  # Track visited nodes
        component_count = 0    # Number of connected components
        for i in range(d):
            if alive[i] and not visited[i]:
                component_count += 1
                bfs = deque([i])
                visited[i] = True
                # Standard BFS to traverse the component
                while bfs:
                    u = bfs.popleft()
                    for w in induced_adj[u]:
                        if alive[w] and not visited[w]:
                            visited[w] = True
                            bfs.append(w)
        # Return the number of connected components in the k-core of the neighborhood
        return component_count

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
        # G: Graph object with vertex_num (number of vertices) and adj_list (adjacency list)
        n = G.vertex_num
        adj = G.adj_list
        sd = [0] * n
        core_cache = dict()  # 缓存：frozenset(neighbors) -> 结果

        for v in range(n):
            neighbors = adj[v]
            if len(neighbors) < k:
                sd[v] = 0
                continue
            key = frozenset(neighbors)
            if key in core_cache:
                sd[v] = core_cache[key]
            else:
                res = kCoreBaseStructuralDiversity._process_single_vertex((adj, v, k))
                core_cache[key] = res
                sd[v] = res
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
