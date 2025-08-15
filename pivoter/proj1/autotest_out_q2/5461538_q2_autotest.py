#!/usr/bin/env python3
# Auto-generated for 5461538

STUDENT_ID = "5461538"
STUDENT_NAME = "Zeyu Zhang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################


class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # compute global core numbers of every vertex with a linear peeling algorithm
    @staticmethod
    def _core_numbers(G) -> list[int]:
        # number of vertices
        num_vertices = G.vertex_num
        # degree list where degree[i] equals the current degree of vertex i
        degree = [len(G.adj_list[v]) for v in range(num_vertices)]
        # result container for core numbers
        core = [0] * num_vertices
        # return early if the graph is empty
        if num_vertices == 0:
            return core
        # counting sort buckets for vertices grouped by degree
        max_degree = max(degree)
        bucket_sizes = [0] * (max_degree + 1)
        for deg in degree:
            bucket_sizes[deg] += 1
        # start index of each degree bucket in the ordered array
        bucket_start = [0] * (max_degree + 1)
        prefix_sum = 0
        for deg in range(max_degree + 1):
            bucket_start[deg] = prefix_sum
            prefix_sum += bucket_sizes[deg]
        # ordered holds vertices sorted by degree from low to high
        ordered = [0] * num_vertices
        # position array stores where each vertex currently sits in ordered
        position = [0] * num_vertices
        for v, deg in enumerate(degree):
            insert_pos = bucket_start[deg]
            ordered[insert_pos] = v
            position[v] = insert_pos
            bucket_start[deg] += 1
        # restore bucket_start so that it again points to bucket starts
        for deg in range(max_degree, 0, -1):
            bucket_start[deg] = bucket_start[deg - 1]
        bucket_start[0] = 0
        # main peeling loop removes vertices from low to high degree
        for i in range(num_vertices):
            v = ordered[i]
            for u in G.adj_list[v]:
                if degree[u] > degree[v]:
                    deg_u = degree[u]
                    index_u = position[u]
                    index_w = bucket_start[deg_u]
                    vertex_w = ordered[index_w]
                    if u != vertex_w:
                        ordered[index_u], ordered[index_w] = vertex_w, u
                        position[u], position[vertex_w] = index_w, index_u
                    bucket_start[deg_u] += 1
                    degree[u] -= 1
            core[v] = degree[v]
        return core

    # breadth first search to extract connected components inside a subset
    @staticmethod
    def _find_components(adj: dict[int, list[int]], vertices: set[int]) -> list[set[int]]:
        # visited set prevents revisiting the same vertex
        visited: set[int] = set()
        # list of all components to return
        comps: list[set[int]] = []
        for v in vertices:
            if v in visited:
                continue
            # start a new component from vertex v
            comp = {v}
            q = deque([v])
            visited.add(v)
            # standard bfs expansion
            while q:
                cur = q.popleft()
                for nbr in adj[cur]:
                    if nbr in vertices and nbr not in visited:
                        visited.add(nbr)
                        q.append(nbr)
                        comp.add(nbr)
            comps.append(comp)
        return comps

    # count connected components inside an arbitrary vertex set
    @staticmethod
    def _count_components(G, vertices: set[int]) -> int:
        visited: set[int] = set()
        components = 0
        for v in vertices:
            if v in visited:
                continue
            components += 1
            stack = [v]
            visited.add(v)
            while stack:
                cur = stack.pop()
                for nbr in G.adj_list[cur]:
                    if nbr in vertices and nbr not in visited:
                        visited.add(nbr)
                        stack.append(nbr)
        return components

    # compute the number of k core connected components inside a vertex set
    @staticmethod
    def _count_k_core_components(G, k, vertices: set[int]) -> int:
        # empty set yields zero components
        if not vertices:
            return 0
        # k equal to zero means every vertex is a zero core
        if k == 0:
            return kCoreBaseStructuralDiversity._count_components(G, vertices)
        # if the set is too small it cannot host a k core
        if len(vertices) < k + 1:
            return 0
        # local degree dictionary restricted to the vertex set
        local_deg = {
            v: sum(1 for u in G.adj_list[v] if u in vertices)
            for v in vertices
        }
        # queue stores vertices that currently violate the k core condition
        queue = deque([v for v, deg in local_deg.items() if deg < k])
        removed: set[int] = set(queue)
        # iterative peeling inside the vertex set
        while queue:
            v = queue.popleft()
            for u in G.adj_list[v]:
                if u in local_deg and u not in removed:
                    local_deg[u] -= 1
                    if local_deg[u] < k:
                        removed.add(u)
                        queue.append(u)
        # remaining vertices form the k core subgraph
        k_core_vertices = vertices - removed
        if not k_core_vertices:
            return 0
        return kCoreBaseStructuralDiversity._count_components(G, k_core_vertices)

    # update tau values for every vertex by iterating over components
    @staticmethod
    def _process_components(G, k, components: list[set[int]], tau: list[int]) -> None:
        # cache avoids recomputing identical neighbor induced subgraphs
        cache: dict[frozenset, int] = {}
        for comp in components:
            # relevant vertices include the component and all its neighbors
            relevant: set[int] = set(comp)
            for v in comp:
                relevant.update(G.adj_list[v])
            for v in relevant:
                nbrs = G.adj_list[v]
                if not nbrs:
                    continue
                # fast path for degree one
                if len(nbrs) == 1:
                    k_core_cc = 1 if k == 0 else 0
                    if k_core_cc > tau[v]:
                        tau[v] = k_core_cc
                    continue
                nbr_set = set(nbrs)
                key = frozenset(nbr_set)
                if key in cache:
                    k_core_cc = cache[key]
                else:
                    k_core_cc = kCoreBaseStructuralDiversity._count_k_core_components(G, k, nbr_set)
                    cache[key] = k_core_cc
                if k_core_cc > tau[v]:
                    tau[v] = k_core_cc

    # main entry that returns tau_k for every vertex
    @staticmethod
    def process(G, k):
        # total number of vertices
        n = G.vertex_num
        # handle empty graph or invalid k
        if n == 0 or k < 0:
            return [0] * n
        # compute global core numbers once
        core_num = kCoreBaseStructuralDiversity._core_numbers(G)
        # filter vertices whose core number is at least k
        valid_vertices = {v for v, c in enumerate(core_num) if c >= k}
        if not valid_vertices:
            return [0] * n
        # extract connected components inside the filtered subgraph
        comps = kCoreBaseStructuralDiversity._find_components(G.adj_list, valid_vertices)
        # initialize output list with zeros
        tau_k = [0] * n
        # fill tau_k by processing each component
        kCoreBaseStructuralDiversity._process_components(G, k, comps, tau_k)
        return tau_k

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
