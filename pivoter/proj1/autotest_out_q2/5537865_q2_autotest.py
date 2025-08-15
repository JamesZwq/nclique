#!/usr/bin/env python3
# Auto-generated for 5537865

STUDENT_ID = "5537865"
STUDENT_NAME = "Jie Gao"

# ======= 学生代码 =======
################################################################################
from collections import deque
from typing import List, Set, Dict

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k: int) -> List[int]:
        n = G.vertex_num
        sd = [0] * n
        if k < 0 or n == 0:
            return sd

        # Step 1: Precompute each vertex’s core number for the whole graph
        core_numbers = kCoreBaseStructuralDiversity._compute_core_numbers(G)

        # Step 2: Build the set of vertices that are in the (k+1)-core
        k1_vertices = {v for v in range(n) if core_numbers[v] >= k+1}
        if not k1_vertices:
            # If there is no (k+1)-core, fall back to a local method for every vertex
            return kCoreBaseStructuralDiversity._fallback_k_core_method(G, k, core_numbers)

        # Step 3: Find all connected components within the (k+1)-core
        k1_components = kCoreBaseStructuralDiversity._find_components(G, k1_vertices)

        # Step 4: Process each (k+1)-core component as a batch
        for comp in k1_components:
            kCoreBaseStructuralDiversity._process_component_batch(G, comp, k, sd)

        # Step 5: For all remaining vertices outside the (k+1)-core, use the local method
        processed = set().union(*k1_components)
        for v in range(n):
            if v not in processed:
                sd[v] = kCoreBaseStructuralDiversity._compute_tau_k_local(
                    G, v, k, core_numbers
                )
        return sd

    @staticmethod
    def _compute_core_numbers(G) -> List[int]:
        """Compute each vertex’s core number in linear time using peeling + counting sort."""
        n = G.vertex_num
        deg = [len(G.adj_list[v]) for v in range(n)]
        maxd = max(deg) if deg else 0

        # Prepare for counting sort by degree
        count = [0] * (maxd + 1)
        for d in deg:
            count[d] += 1
        start = [0] * (maxd + 1)
        for i in range(1, maxd + 1):
            start[i] = start[i-1] + count[i-1]

        # Sort vertices by degree
        pos = list(start)
        vert = [0] * n
        for v in range(n):
            d = deg[v]
            vert[pos[d]] = v
            pos[d] += 1

        # Perform the peeling process
        core_number = [0] * n
        for i in range(n):
            u = vert[i]
            core_number[u] = deg[u]
            for w in G.adj_list[u]:
                if deg[w] > deg[u]:
                    deg[w] -= 1
                    # We could update 'pos' here for exact ordering, but it's not required for correctness
        return core_number

    @staticmethod
    def _find_components(G, vertices: Set[int]) -> List[List[int]]:
        """Use BFS to find all connected components within the given subset of vertices."""
        visited = set()
        comps = []
        for v in vertices:
            if v in visited:
                continue
            queue = deque([v])
            visited.add(v)
            comp = []
            while queue:
                u = queue.popleft()
                comp.append(u)
                for w in G.adj_list[u]:
                    if w in vertices and w not in visited:
                        visited.add(w)
                        queue.append(w)
            comps.append(comp)
        return comps

    @staticmethod
    def _process_component_batch(G, component: List[int], k: int, sd: List[int]):
        """
        For a single (k+1)-core component, compute τ_k for each neighbor v
        based on k-core counts within that component.
        """
        comp_set = set(component)
        # Collect all neighbors of the component vertices
        all_nbrs = set()
        for u in component:
            all_nbrs.update(G.adj_list[u])

        for v in all_nbrs:
            nbrs = set(G.adj_list[v]) & comp_set
            if len(nbrs) < k:
                continue
            # Use the fast local k-core count on this neighbor set
            sd[v] += kCoreBaseStructuralDiversity._fast_k_core_count(
                G, nbrs, k, comp_set
            )

    @staticmethod
    def _fast_k_core_count(G, candidates: Set[int], k: int, valid_neigh: Set[int]) -> int:
        """
        Quickly peel vertices of degree < k from the induced subgraph,
        then count how many connected sub-cores remain.
        """
        vertices = set(candidates)
        # 1) Remove all vertices with current degree < k
        deg = {}
        for u in vertices:
            deg[u] = sum(1 for w in G.adj_list[u] if w in vertices and w in valid_neigh)
        queue = deque(u for u in vertices if deg[u] < k)
        removed = set(queue)
        while queue:
            u = queue.popleft()
            for w in G.adj_list[u]:
                if w in vertices and w not in removed and w in valid_neigh:
                    deg[w] -= 1
                    if deg[w] < k:
                        removed.add(w)
                        queue.append(w)
        core_vs = vertices - removed
        if not core_vs:
            return 0

        # 2) Count connected components in the remaining subgraph
        visited = set()
        cnt = 0
        for u in core_vs:
            if u in visited:
                continue
            cnt += 1
            dq = deque([u])
            visited.add(u)
            while dq:
                x = dq.popleft()
                for w in G.adj_list[x]:
                    if w in core_vs and w not in visited:
                        visited.add(w)
                        dq.append(w)
        return cnt

    @staticmethod
    def _compute_tau_k_local(G, v: int, k: int, core_numbers: List[int]) -> int:
        """
        Fallback method: for a single vertex v,
        build its neighbor-induced subgraph (neighbors with core_number ≥ k)
        and compute τ_k(v) exactly.
        """
        nbrs = [u for u in G.adj_list[v] if core_numbers[u] >= k]
        if len(nbrs) < k:
            return 0

        # Build adjacency for the induced subgraph
        subg = {u: [] for u in nbrs}
        for u in nbrs:
            for w in G.adj_list[u]:
                if w in subg:
                    subg[u].append(w)

        # Reuse fast k-core count on this subgraph
        return kCoreBaseStructuralDiversity._fast_k_core_count(
            G, set(nbrs), k, set(nbrs)
        )

    @staticmethod
    def _fallback_k_core_method(G, k: int, core_numbers: List[int]) -> List[int]:
        """
        If there is no (k+1)-core, fall back to computing τ_k(v)
        locally for every vertex v in the graph.
        """
        n = G.vertex_num
        sd = [0] * n
        for v in range(n):
            sd[v] = kCoreBaseStructuralDiversity._compute_tau_k_local(
                G, v, k, core_numbers
            )
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
