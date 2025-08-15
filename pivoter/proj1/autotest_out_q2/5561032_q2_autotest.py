#!/usr/bin/env python3
# Auto-generated for 5561032

STUDENT_ID = "5561032"
STUDENT_NAME = "Chiyue Rao"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k): # O(V*E)
        n = G.vertex_num
        # compute the core numbers of all vertices
        core = kCoreBaseStructuralDiversity._compute_core_numbers(G.adj_list, n)

        sd = [0] * n
        for v in range(n): # O(V)
            neighbors = G.adj_list[v]
            if not neighbors:
                continue

            # only consider neighbors with core number >= k, which is faster than checking all neighbors
            nbr_set = {u for u in neighbors if core[u] >= k} # O(E)
            if not nbr_set:
                continue

            # build the neighbor-induced subgraph (only contains nbr_set), all the neighbors' degree are >= k
            nbr_adj = {u: [w for w in G.adj_list[u] if w in nbr_set] for u in nbr_set} # O(E)

            # find the connected components in the k-core
            kcore_nodes = kCoreBaseStructuralDiversity._compute_k_core(nbr_adj, k) # O(E)
            components = kCoreBaseStructuralDiversity._find_connected_components(nbr_adj, kcore_nodes) # O(E)
            sd[v] = len(components)

        return sd

    # Flat Array k-core algorithm
    @staticmethod
    def _compute_core_numbers(adj_list, n): # O(V + E)
        degree = [len(adj_list[u]) for u in range(n)] # O(V)
        max_deg = max(degree) if n > 0 else 0

        # bin[d] = a list with the number of vertices with degree d (d is the index)
        bin = [0] * (max_deg + 1)
        for d in degree: # O(V)
            bin[d] += 1

        # start[d] = the starting position of vertices with degree d in the sorted array
        start = [0] * (max_deg + 1)
        for d in range(1, max_deg + 1): # O(V)
            # start[d-1] is the starting position of vertices with degree d-1 in the sorted array
            # bin[d-1] is the number of vertices with degree d-1
            start[d] = start[d - 1] + bin[d - 1] 

        # vert[] = the sorted array of vertices by degree
        vert = [0] * n
        pos = [0] * n
        for v in range(n): # O(V)
            d = degree[v]
            # pos[v] is the starting position of vertex v in the sorted array
            pos[v] = start[d]
            # vert[start[d]] = v
            vert[pos[v]] = v
            start[d] += 1

        # correct start
        for d in range(max_deg, 0, -1): # O(V)
            start[d] = start[d - 1]
        start[0] = 0

        # main loop: update the core number, the key trick of the flat array k-core algorithm
        # For each neighbor u of vertex v, if u has a higher degree than v,
        # we decrease u's degree because v is being removed.
        # To maintain vertices sorted by their degrees without re-sorting,
        # we use three arrays:
        #   pos[u]   -> current position of u in the degree-sorted array
        #   vert[i]  -> vertex currently at position i
        #   start[d] -> starting index of vertices with degree d
        # When u's degree changes, we swap u with the vertex w at the start
        # position of its degree bucket and update their positions.
        # This ensures that vertices remain grouped by degree efficiently,
        # which is the key trick of the flat array k-core algorithm.
        for i in range(n): # O(V)
            v = vert[i]
            for u in adj_list[v]: # O(E)
                if degree[u] > degree[v]:
                    du = degree[u]
                    pu = pos[u]
                    pw = start[du]
                    w = vert[pw]
                    if u != w:
                        pos[u], pos[w] = pw, pu
                        vert[pu], vert[pw] = w, u
                    start[du] += 1
                    degree[u] -= 1
        return degree

    # k-core peeling algorithm (for the subgraph)
    @staticmethod
    def _compute_k_core(adj, k): # O(E)
        degree = {u: len(adj[u]) for u in adj} # O(E)
        queue = deque([u for u in adj if degree[u] < k]) # O(E)
        removed = set()
        while queue: # O(E)
            u = queue.popleft()
            removed.add(u)
            for v in adj[u]: # O(E)
                if v not in removed:
                    degree[v] -= 1
                    if degree[v] == k - 1:
                        queue.append(v)
        return set(adj.keys()) - removed

    # connected components
    @staticmethod
    def _find_connected_components(adj, vertices): # O(E)
        visited = set()
        components = []
        for u in vertices: # O(V)
            if u in visited:
                continue
            queue = deque([u])
            visited.add(u)
            component = [u]
            while queue: # O(E)
                cur = queue.popleft()
                for nei in adj[cur]: # O(E)
                    if nei in vertices and nei not in visited:
                        visited.add(nei)
                        queue.append(nei)
                        component.append(nei)
            components.append(component)
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
