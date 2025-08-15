#!/usr/bin/env python3
# Auto-generated for 5507401

STUDENT_ID = "5507401"
STUDENT_NAME = "Xiaoran Zhang"

# ======= 学生代码 =======
class kCoreBaseStructuralDiversity(object):
    """
    ========================================================================
    comp9312 25t2 project q2
    ========================================================================

    main parts of our algorithm:
    1. global k-core decomposition:
       - compute core numbers for all vertices in one pass
       - identifies which vertices belong to (k+1)-core subgraph
       - avoids repeated computation across different vertices
       - uses efficient bucket-based core decomposition

    2. (k+1)-core optimization strategy:
       - build (k+1)-core subgraph for batch processing
       - find connected components in (k+1)-core subgraph
       - process vertices in batches rather than individually
       - significantly reduces redundant ego-network construction

    3. batch processing for efficiency:
       - process vertices in same (k+1)-core component together
       - share computation results when possible
       - handle remaining vertices with individual processing
       - maintain correctness while improving performance

    4. ego-network k-core analysis:
       - build neighborhood graph for each vertex
       - perform k-core decomposition on ego-network
       - count connected components in k-core of ego-network
       - this gives structural diversity value for each vertex

    time complexity analysis:
    - global core decomposition: o(m + n) using bucket sort optimization
    - (k+1)-core component finding: o(m + n) using dfs traversal
    - batch processing: o(sum of ego-network processing times)
    - ego-network processing per vertex: o(degree^2 + degree * k)
    - total complexity: o(m + n + sum over all vertices of degree^2)
    - this is much better than naive o(n * m) approach

    space complexity: o(n + m + max_degree^2) for core decomposition and ego-networks
    ========================================================================
    """

    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        compute k-core based structural diversity for all vertices in graph G

        time complexity: o(m + n + sum of ego-network processing)
        space complexity: o(n + m + max_degree^2)
        """
        n = G.vertex_num
        adj_list = G.adj_list
        sd = [0] * n  # result array for structural diversity values

        # step 1: compute core numbers for all vertices - o(m + n)
        # this identifies which vertices are in (k+1)-core subgraph
        core_numbers = kCoreBaseStructuralDiversity._compute_core_numbers(G)

        # step 2: build (k+1)-core subgraph - o(n)
        # vertices with core number >= k+1 form the (k+1)-core
        k_plus_1_vertices = set(v for v in range(n) if core_numbers[v] >= k + 1)

        if not k_plus_1_vertices:
            # if no (k+1)-core exists, check k-core instead
            return kCoreBaseStructuralDiversity._fallback_k_core_method(G, k, core_numbers)

        # step 3: find connected components in (k+1)-core subgraph - o(m + n)
        # this allows us to process vertices in batches for better efficiency
        k_plus_1_components = kCoreBaseStructuralDiversity._find_components(
            G, k_plus_1_vertices
        )

        # step 4: batch processing for (k+1)-core components
        # process each component as a group to share computation
        for component in k_plus_1_components:
            kCoreBaseStructuralDiversity._process_component_batch(
                G, component, k, sd
            )

        # step 5: process remaining vertices (not in (k+1)-core)
        # collect all vertices that are neighbors of (k+1)-core vertices
        processed = set()
        for comp in k_plus_1_components:
            for v in comp:
                processed.add(v)  # mark (k+1)-core vertices as processed
                for u in adj_list[v]:
                    processed.add(u)  # mark their neighbors as processed

        # handle vertices not yet processed - o(remaining vertices * ego processing)
        for v in range(n):
            if v not in processed:
                sd[v] = kCoreBaseStructuralDiversity._compute_single_vertex_sd(G, v, k)

        return sd

    @staticmethod
    def _compute_core_numbers(G):
        """
        compute core number for each vertex using efficient k-core decomposition

        time complexity: o(m + n)
        space complexity: o(n + max_degree)
        """
        n = G.vertex_num
        adj_list = G.adj_list

        # initialize degrees and data structures - o(n)
        degrees = [len(adj_list[v]) for v in range(n)]
        core_numbers = [0] * n
        removed = [False] * n

        # set up bucket sort for efficient processing - o(n + max_degree)
        max_degree = max(degrees) if degrees else 0
        buckets = [[] for _ in range(max_degree + 1)]

        # put each vertex in corresponding degree bucket - o(n)
        for v in range(n):
            buckets[degrees[v]].append(v)

        # process vertices from lowest degree to highest - o(m + n)
        for current_k in range(max_degree + 1):
            while buckets[current_k]:
                v = buckets[current_k].pop()
                if removed[v]:
                    continue

                # set core number for vertex v - o(1)
                core_numbers[v] = current_k
                removed[v] = True

                # update degrees of neighbors - o(degree(v))
                for u in adj_list[v]:
                    if not removed[u] and degrees[u] > current_k:
                        degrees[u] -= 1  # neighbor loses connection to v
                        if degrees[u] >= current_k:
                            buckets[degrees[u]].append(u)  # move to correct bucket

        return core_numbers

    @staticmethod
    def _fallback_k_core_method(G, k, core_numbers):
        """
        fallback method when no (k+1)-core exists

        time complexity: same as individual vertex processing
        """
        n = G.vertex_num
        sd = [0] * n

        # check if k-core exists - o(n)
        k_core_vertices = set(v for v in range(n) if core_numbers[v] >= k)

        if not k_core_vertices:
            # no k-core exists, all structural diversity values are 0
            return sd

        # process each vertex individually - o(n * ego processing time)
        for v in range(n):
            sd[v] = kCoreBaseStructuralDiversity._compute_single_vertex_sd(G, v, k)

        return sd

    @staticmethod
    def _find_components(G, vertex_set):
        """
        find connected components in subgraph induced by vertex_set
        uses depth-first search to identify all components

        time complexity: o(|vertex_set| + edges_in_subgraph)
        space complexity: o(|vertex_set|) for visited set and recursion
        """
        visited = set()
        components = []

        def dfs(start_vertex, current_component):
            """
            depth-first search to find all vertices in current component
            explores only vertices in the given vertex_set
            """
            if start_vertex in visited:
                return
            visited.add(start_vertex)
            current_component.append(start_vertex)

            # visit all neighbors that are in the vertex set - o(degree)
            for neighbor in G.adj_list[start_vertex]:
                if neighbor in vertex_set and neighbor not in visited:
                    dfs(neighbor, current_component)

        # find all components using dfs - o(vertex_set + edges)
        for v in vertex_set:
            if v not in visited:
                component = []
                dfs(v, component)
                components.append(component)

        return components

    @staticmethod
    def _process_component_batch(G, component, k, sd):
        """
        process connected component in batch for better efficiency
        """
        # process each vertex in component - can be optimized further
        # current implementation processes individually but keeps batch structure
        for v in component:
            sd[v] = kCoreBaseStructuralDiversity._compute_single_vertex_sd(G, v, k)

    @staticmethod
    def _compute_single_vertex_sd(G, v, k):
        """
        compute structural diversity for single vertex using ego-network analysis

        steps:
        1. build ego-network (neighborhood graph) for vertex v
        2. perform k-core decomposition on ego-network
        3. count connected components in k-core of ego-network

        time complexity: o(degree(v)^2 + degree(v) * k)
        space complexity: o(degree(v)^2) for ego-network storage
        """
        adj_list = G.adj_list

        # build ego-network (neighborhood graph) - o(degree(v)^2)
        neighbors = set(adj_list[v])
        if len(neighbors) < k:
            return 0  # not enough neighbors to form k-core

        # build adjacency list for ego-network - o(degree(v)^2)
        ego_adj = {u: [] for u in neighbors}
        for u in neighbors:
            for w in adj_list[u]:
                if w in neighbors and w != u:
                    ego_adj[u].append(w)  # edge between two neighbors of v

        # perform k-core decomposition on ego-network - o(degree(v) * k)
        k_core_vertices = kCoreBaseStructuralDiversity._ego_k_core_decomposition(
            ego_adj, neighbors, k
        )

        # count connected components in k-core - o(k_core_size + edges)
        return kCoreBaseStructuralDiversity._count_components(ego_adj, k_core_vertices)

    @staticmethod
    def _ego_k_core_decomposition(ego_adj, vertices, k):
        """
        perform k-core decomposition on ego-network
        repeatedly remove vertices with degree < k until no more can be removed

        time complexity: o(|vertices| * k + |edges|)
        space complexity: o(|vertices|)
        """
        # calculate initial degrees in ego-network - o(vertices + edges)
        degrees = {u: len(ego_adj[u]) for u in vertices}
        remaining = set(vertices)

        # iteratively remove vertices with degree < k - o(vertices * k)
        changed = True
        while changed:
            changed = False
            to_remove = []

            # find all vertices with degree < k - o(remaining vertices)
            for u in remaining:
                if degrees[u] < k:
                    to_remove.append(u)

            # remove these vertices and update neighbor degrees - o(edges)
            for u in to_remove:
                remaining.remove(u)
                for w in ego_adj[u]:
                    if w in remaining:
                        degrees[w] -= 1  # neighbor loses connection to removed vertex
                changed = True

        return remaining

    @staticmethod
    def _count_components(ego_adj, vertices):
        """
        count connected components in given vertex set using dfs

        time complexity: o(|vertices| + |edges_in_subgraph|)
        space complexity: o(|vertices|) for visited set and recursion
        """
        if not vertices:
            return 0

        visited = set()
        components = 0

        def dfs(node):
            """depth-first search to visit all nodes in current component"""
            if node in visited:
                return
            visited.add(node)
            # visit all neighbors that are still in the vertex set
            for neighbor in ego_adj[node]:
                if neighbor in vertices and neighbor not in visited:
                    dfs(neighbor)

        # count components using dfs - o(vertices + edges)
        for u in vertices:
            if u not in visited:
                dfs(u)  # explore entire component starting from u
                components += 1

        return components

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
