#!/usr/bin/env python3
# Auto-generated for 5556885

STUDENT_ID = "5556885"
STUDENT_NAME = "Thomas  Jian Shao"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def compute_coreness(G):
        """Standard k-core decomposition: returns the coreness value for each vertex."""
        
        # n = number of vertices in the graph
        n = G.vertex_num
        # deg[v] = degree of vertex v (number of neighbors)
        deg = [len(G.adj_list[v]) for v in range(n)]
        # coreness[v] will store the coreness value of vertex v
        coreness = [0] * n
        # Bucket sort structure: bucket[d] stores vertices with current degree d
        # The maximum possible degree is max(deg), so we need that many buckets
        bucket = [[] for _ in range(max(deg) + 1)]
    
        # Initialize the buckets by placing each vertex into the bucket of its degree
        for v in range(n):
            bucket[deg[v]].append(v)
    
        # curr_k tracks the current core level we are peeling at
        curr_k = 0
        # removed counts how many vertices have been processed
        removed = 0
    
        # Process vertices in increasing degree order
        for d in range(len(bucket)):
            # While there are vertices in the current degree bucket
            while bucket[d]:
                # Remove one vertex with current degree d
                v = bucket[d].pop()
                # If the degree of v has changed since it was bucketed, skip it
                if deg[v] != d:
                    continue
                # Assign the coreness of v as its current degree d
                coreness[v] = d
                removed += 1
    
                # For each neighbor u of v
                for u in G.adj_list[v]:
                    # If u has degree greater than d, we need to decrement its degree
                    # because v is being removed
                    if deg[u] > d:
                        deg[u] -= 1
                        # Place u into the new bucket corresponding to its updated degree
                        bucket[deg[u]].append(u)
    
        # Return the coreness values for all vertices
        return coreness


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
        sd = [0] * n

        # Step 0: Global coreness decomposition
        # Compute the global coreness of each vertex using k-core decomposition.
        # This value indicates the highest k for which the vertex belongs to a k-core globally.
        coreness = kCoreBaseStructuralDiversity.compute_coreness(G)

        for v in range(n):
            # Step 1: Candidate neighbor set
            # Select neighbors of v whose global coreness is >= k.
            # These are the only candidates that could possibly survive in the k-core
            # of the neighbor-induced subgraph.
            neighbors = set(u for u in G.adj_list[v] if coreness[u] >= k)
            if not neighbors:
                # If there are no valid neighbors, the diversity is zero.
                continue

            # Step 2: Initialize degree counts in the neighbor-induced subgraph
            # For each candidate neighbor, count how many connections it has 
            # within the candidate neighbor set.
            degree = {u: 0 for u in neighbors}
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbors:
                        degree[u] += 1

            # Step 3: Peeling process (k-core pruning in neighbor subgraph)
            # Remove nodes with degree < k iteratively, adjusting the degree
            # of their neighbors, until all remaining nodes satisfy degree >= k.
            queue = deque([u for u in neighbors if degree[u] < k])
            while queue:
                u = queue.popleft()
                for w in G.adj_list[u]:
                    if w in degree and degree[w] >= k:
                        degree[w] -= 1
                        # If w’s degree falls below k after removing u,
                        # enqueue it for removal as well.
                        if degree[w] == k-1:
                            queue.append(w)
                degree[u] = -1  # Mark u as removed

            # Step 4: Count connected components in the remaining k-core
            # After pruning, the vertices with degree >= k form the k-core.
            # We now count how many connected components exist among them.
            visited = set()
            count = 0
            for u in neighbors:
                if degree[u] >= k and u not in visited:
                    # Found a new component
                    count += 1
                    stack = [u]
                    visited.add(u)
                    # Depth-first search to mark all nodes in this component
                    while stack:
                        cur = stack.pop()
                        for w in G.adj_list[cur]:
                            if w in degree and degree[w] >= k and w not in visited:
                                visited.add(w)
                                stack.append(w)
            # τ_k(v) = number of k-core connected components in v’s neighbor-induced subgraph
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
