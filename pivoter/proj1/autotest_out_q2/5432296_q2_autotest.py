#!/usr/bin/env python3
# Auto-generated for 5432296

STUDENT_ID = "5432296"
STUDENT_NAME = "Zihao Yan"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    """
    Provides methods to compute the k-core based structural diversity metric for all vertices in a graph.
    This metric, tau_k(v), for a vertex v is defined as the number of connected components
    in the k-core of the subgraph induced by the neighbors of v.
    """
    def __init__(self):
        pass

    @staticmethod
    def compute_core_numbers(G):
        """
        Computes the core number for each vertex in the graph G using a linear-time
        core decomposition algorithm. This algorithm has a time complexity of O(n + m),
        where n is the number of vertices and m is the number of edges.

        Returns:
            List[int]: A list where the i-th element is the core number of vertex i.
        """
        n = G.vertex_num
        # Store the initial degrees of all vertices.
        degree = [len(G.adj_list[v]) for v in range(n)]
        max_d = max(degree) if n else 0

        # --- Degree Sorting Initialization ---
        # The following block sorts vertices by their degree in linear time,
        # similar to counting sort.
        # `bin_count[d]` will store the number of vertices with degree d.
        bin_count = [0] * (max_d + 1)
        for d in degree:
            bin_count[d] += 1
        
        # `start[d]` will store the starting position for vertices of degree d in the sorted array.
        start = [0] * (max_d + 1)
        acc = 0
        for d in range(max_d + 1):
            start[d], acc = acc, acc + bin_count[d]

        # `vert` is the array of vertices sorted by their current degree.
        # `pos[v]` stores the position of vertex v in the `vert` array.
        vert = [0] * n
        pos = [0] * n
        for v in range(n):
            d = degree[v]
            vert[start[d]] = v
            pos[v] = start[d]
            start[d] += 1

        # Restore `start` array to point to the beginning of each degree bin.
        for d in range(max_d, 0, -1):
            start[d] = start[d-1]
        start[0] = 0

        cores = [0] * n
        # --- Peeling Process ---
        # Process vertices in increasing order of their degree.
        for i in range(n):
            v = vert[i]
            cores[v] = degree[v]
            
            # For each neighbor u of v, if its degree is higher, decrement it.
            for u in G.adj_list[v]:
                if degree[u] > degree[v]:
                    du = degree[u]
                    pu = pos[u]
                    
                    # To maintain the sorted order in `vert`, swap u with the first vertex
                    # `w` in the bin for degree `du`.
                    pw = start[du]
                    w = vert[pw]
                    
                    # Perform the swap in `vert` and update positions in `pos`.
                    vert[pu], vert[pw] = w, u
                    pos[w], pos[u] = pu, pw
                    
                    # Increment the start of the bin for degree `du` as it now contains one less vertex.
                    start[du] += 1
                    # Decrement the degree of u.
                    degree[u] -= 1
        return cores

    @staticmethod
    def compute_k_core_count(V_set, G, k):
        """
        Builds a subgraph from a given vertex set, peels it to find the k-core,
        and returns the number of connected components in the resulting k-core.
        """
        m = len(V_set)
        if m <= k:  # A k-core requires at least k+1 vertices, but check for m<=k is safer.
            return 0

        # Map original vertex IDs to local indices (0 to m-1) for efficient processing.
        idx = {v:i for i, v in enumerate(V_set)}
        # Build the adjacency list and degree array for the subgraph.
        adj = [[] for _ in range(m)]
        deg = [0] * m
        for i, v in enumerate(V_set):
            for u in G.adj_list[v]:
                if u in idx:
                    adj[i].append(idx[u])
            deg[i] = len(adj[i])

        # Initialize a queue with vertices whose degree is less than k for peeling.
        queue = deque()
        removed = [False] * m
        for i, d in enumerate(deg):
            if d < k:
                removed[i] = True
                queue.append(i)

        # Peeling process: iteratively remove vertices with degree < k.
        while queue:
            x = queue.popleft()
            for y in adj[x]:
                if not removed[y]:
                    deg[y] -= 1
                    if deg[y] < k:
                        removed[y] = True
                        queue.append(y)

        # Count connected components in the remaining graph (the k-core).
        visited = [False] * m
        cnt = 0
        for i in range(m):
            if not removed[i] and not visited[i]:
                cnt += 1
                dq = deque([i])
                visited[i] = True
                # Use BFS to find all vertices in the current component.
                while dq:
                    cur = dq.popleft()
                    for nb in adj[cur]:
                        if not removed[nb] and not visited[nb]:
                            visited[nb] = True
                            dq.append(nb)
        return cnt

    @staticmethod
    def process(G, k):
        """
        Calculates the k-core based structural diversity metric (tau_k) for all vertices in graph G.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The input graph object. It is expected to have attributes `vertex_num` and `adj_list`.
        k : int
            The parameter 'k' for the k-core definition.

        Returns
        -------
        List[int]
            A list where the value at index `v` is the structural diversity metric tau_k(v).
        """
        n = G.vertex_num
        if k < 0 or n == 0:
            return [0] * n
        # If k is larger than the number of vertices, no k-core can exist.
        if k > n:
            return [0] * n

        # 1. Perform global k-core decomposition to get core numbers for all vertices.
        cores = kCoreBaseStructuralDiversity.compute_core_numbers(G)
        tau = [0] * n

        # 2. Identify candidate vertices for computation (those in the (k+1)-core).
        candidates = [v for v in range(n) if cores[v] >= k+1]
        if not candidates:
            return tau

        # Build a subgraph containing only the candidate vertices and edges between them.
        sub_adj = {v: [u for u in G.adj_list[v] if cores[u] >= k+1]
                   for v in candidates}

        # 3. Find connected components within the candidate subgraph.
        visited = set()
        for v in candidates:
            if v in visited:
                continue
            
            # Use BFS to find all vertices in one connected component `C`.
            C = []
            dq = deque([v])
            visited.add(v)
            while dq:
                x = dq.popleft()
                C.append(x)
                for u in sub_adj[x]:
                    if u not in visited:
                        visited.add(u)
                        dq.append(u)

            # 4. For the component C, define a `related` set of vertices to compute tau_k for.
            #    related = (vertices in C) U (neighbors of vertices in C).
            related = set(C)
            for x in C:
                related.update(G.adj_list[x])

            # For each vertex `u` in the `related` set, compute its tau_k value.
            for u in related:
                # The subgraph is induced by neighbors of u with a core number of at least k.
                nbrs = [w for w in G.adj_list[u] if cores[w] >= k]
                cnt = kCoreBaseStructuralDiversity.compute_k_core_count(nbrs, G, k)
                # Update tau[u] if a larger component count is found.
                if cnt > tau[u]:
                    tau[u] = cnt

        return tau


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    #
    # Algorithm Description and Complexity Analysis:
    #
    # The algorithm computes tau_k(v), the number of connected components in the
    # k-core of the subgraph induced by v's neighbors.
    #
    # 1. Global Core Decomposition: First, it computes the core number for all
    #    vertices in O(n+m) using a linear-time peeling algorithm.
    #
    # 2. Candidate Identification: It identifies connected components of vertices
    #    with a core number of at least k+1. This is the main pool of vertices
    #    that require computation.
    #
    # 3. Local Diversity Calculation: For each component and its neighbors, it
    #    iterates through each vertex `u`. For each `u`, it constructs the
    #    subgraph induced by its neighbors and calls `compute_k_core_count`
    #    to find the number of components in its local k-core.
    #
    # Complexity:
    # - The pre-computation (global core decomposition, candidate finding) takes O(n+m).
    # - The main bottleneck is the final calculation loop. For each processed
    #   vertex `u`, building its local subgraph takes time proportional to the sum of
    #   the degrees of its neighbors.
    # - This leads to an overall complexity that is bounded by O(n*Delta^2 + m),
    #   where n is the number of vertices, m is the number of edges, and Delta
    #   is the maximum degree of the graph.
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
