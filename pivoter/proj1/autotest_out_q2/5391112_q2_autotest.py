#!/usr/bin/env python3
# Auto-generated for 5391112

STUDENT_ID = "5391112"
STUDENT_NAME = "Yong Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from typing import List, Set, Dict
################################################################################


class kCoreBaseStructuralDiversity(object):
    """
    Compute τ_k(v) — k-core-based structural diversity
    Works with UndirectedUnweightedGraph.adj_list only.
    """

    def __init__(self):
        pass

    # -------------------------------------------------------------------------
    # Public entry point
    # -------------------------------------------------------------------------
    @staticmethod
    def process(G, k: int) -> List[int]:
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph   # simple undirected graph
        k : int                        # k-core threshold

        Returns
        -------
        List[int]                      # τ_k(v) for each vertex v
        """
        num_vertices = G.vertex_num
        tau = [0] * num_vertices

        # Trivial cases
        if k <= 0 or num_vertices == 0:
            return tau

        # 1) Core number for every vertex (linear-time peeling)
        core_num = kCoreBaseStructuralDiversity._compute_core_numbers(G)

        # 2) Vertices whose core number ≥ k+1
        upper_core_vertices = {v for v in range(num_vertices)
                               if core_num[v] >= k + 1}
        if not upper_core_vertices:
            # No (k+1)-core exists; fall back to pure k-core computation
            return kCoreBaseStructuralDiversity._fallback_k_core_method(
                G, k, core_num
            )

        # 3) Connected components inside the (k+1)-core subgraph
        upper_components = kCoreBaseStructuralDiversity._find_components(
            G, upper_core_vertices
        )

        # 4) Batch process every (k+1)-core component
        for comp in upper_components:
            kCoreBaseStructuralDiversity._process_component_batch(
                G, comp, k, tau
            )

        # 5) Local computation for vertices not yet covered
        covered: Set[int] = set()
        for comp in upper_components:
            for v in comp:
                covered.add(v)
                covered.update(G.adj_list[v])        # neighbors also considered

        for v in range(num_vertices):
            if v not in covered or tau[v] == 0:
                tau[v] += kCoreBaseStructuralDiversity._compute_tau_k_local(
                    G, v, k, core_num
                )

        return tau

    # -------------------------------------------------------------------------
    # Core number computation
    # -------------------------------------------------------------------------
    @staticmethod
    def _compute_core_numbers(G) -> List[int]:
        """
        Linear-time k-core peeling using bucket sorting.
        Returns a list where entry i is the core number of vertex i.
        """
        n = G.vertex_num
        if n == 0:
            return []

        deg = [len(G.adj_list[v]) for v in range(n)]
        max_deg = max(deg)
        bucket_start = [0] * (max_deg + 1)

        # Histogram of degrees
        for d in deg:
            bucket_start[d] += 1

        # Prefix sum: bucket_start[d] → start index for degree d
        start_idx = 0
        for d in range(max_deg + 1):
            cnt = bucket_start[d]
            bucket_start[d] = start_idx
            start_idx += cnt

        # Bucketed vertex arrays
        ordered = [0] * n
        pos = [0] * n
        for v in range(n):
            d = deg[v]
            pos[v] = bucket_start[d]
            ordered[pos[v]] = v
            bucket_start[d] += 1

        # Restore bucket starts
        for d in range(max_deg, 0, -1):
            bucket_start[d] = bucket_start[d - 1]
        bucket_start[0] = 0

        # Peeling
        for i in range(n):
            v = ordered[i]
            for u in G.adj_list[v]:
                if deg[u] > deg[v]:
                    du = deg[u]
                    pu = pos[u]
                    pw = bucket_start[du]
                    w = ordered[pw]

                    if u != w:
                        pos[u], pos[w] = pw, pu
                        ordered[pu], ordered[pw] = w, u
                    bucket_start[du] += 1
                    deg[u] -= 1

        return deg  # final degrees are core numbers

    # -------------------------------------------------------------------------
    # Connected components inside a vertex subset
    # -------------------------------------------------------------------------
    @staticmethod
    def _find_components(G, vertices: Set[int]) -> List[List[int]]:
        """BFS to enumerate connected components within 'vertices'."""
        seen: Set[int] = set()
        comps: List[List[int]] = []

        for v in vertices:
            if v not in seen:
                comp = []
                q = deque([v])
                seen.add(v)
                while q:
                    u = q.popleft()
                    comp.append(u)
                    for w in G.adj_list[u]:
                        if w in vertices and w not in seen:
                            seen.add(w)
                            q.append(w)
                comps.append(comp)
        return comps

    # -------------------------------------------------------------------------
    # Helpers for (k+1)-core component processing
    # -------------------------------------------------------------------------
    @staticmethod
    def _get_candidate_centers(G, comp_set: Set[int]) -> Set[int]:
        """
        Return all vertices that are 1-hop neighbors of any vertex in comp_set.
        Only these vertices can have ≥ k neighbors inside this component.
        """
        centers: Set[int] = set()
        for v in comp_set:
            centers.update(G.adj_list[v])
        return centers

    @staticmethod
    def _center_contribution(G, center: int,
                             comp_set: Set[int], k: int) -> int:
        """
        How many k-core connected components are present in the subgraph
        induced by 'center's neighbors that lie in comp_set?
        """
        neigh_in_comp = {u for u in G.adj_list[center] if u in comp_set}
        if len(neigh_in_comp) < k:
            return 0
        return kCoreBaseStructuralDiversity._fast_k_core_count(
            G, neigh_in_comp, k
        )

    @staticmethod
    def _process_component_batch(G, component: List[int],
                                 k: int, tau: List[int]) -> None:
        """
        Update tau[] for every valid center vertex related to one
        (k+1)-core connected component.
        """
        comp_set = set(component)
        centers = kCoreBaseStructuralDiversity._get_candidate_centers(
            G, comp_set
        )

        for c in centers:
            tau[c] += kCoreBaseStructuralDiversity._center_contribution(
                G, c, comp_set, k
            )

    # -------------------------------------------------------------------------
    # Fast k-core counting in an induced subgraph
    # -------------------------------------------------------------------------
    @staticmethod
    def _peel_vertices(G, vertices: Set[int], k: int) -> Set[int]:
        """
        Perform k-core peeling on the subgraph induced by 'vertices'.
        Return the surviving vertices (each has degree ≥ k within the subset).
        """
        deg_map: Dict[int, int] = {
            v: sum(1 for u in G.adj_list[v] if u in vertices)
            for v in vertices
        }

        queue = deque([v for v, d in deg_map.items() if d < k])
        removed: Set[int] = set(queue)

        while queue:
            v = queue.popleft()
            for u in G.adj_list[v]:
                if u in vertices and u not in removed:
                    deg_map[u] -= 1
                    if deg_map[u] < k:
                        removed.add(u)
                        queue.append(u)

        return vertices - removed

    @staticmethod
    def _fast_k_core_count(G, vertices: Set[int], k: int) -> int:
        """
        Quickly peel vertices with degree < k, then count components.
        """
        if len(vertices) < k:
            return 0

        remaining = kCoreBaseStructuralDiversity._peel_vertices(
            G, vertices, k
        )
        if not remaining:
            return 0

        return kCoreBaseStructuralDiversity._count_components_fast(
            G, remaining
        )

    # -------------------------------------------------------------------------
    # Connected component counting helpers
    # -------------------------------------------------------------------------
    @staticmethod
    def _count_components_fast(G, vertices: Set[int]) -> int:
        """DFS to count connected components in an induced subgraph."""
        if len(vertices) <= 1:
            return len(vertices)

        seen: Set[int] = set()
        cnt = 0
        for start in list(vertices):
            if start not in seen:
                stack = [start]
                seen.add(start)
                while stack:
                    u = stack.pop()
                    for w in G.adj_list[u]:
                        if w in vertices and w not in seen:
                            stack.append(w)
                            seen.add(w)
                cnt += 1
        return cnt

    # -------------------------------------------------------------------------
    # Local τ_k(v) when batch processing did not cover the vertex
    # -------------------------------------------------------------------------
    @staticmethod
    def _compute_tau_k_local(G, v: int, k: int,
                             core_num: List[int]) -> int:
        """Compute τ_k(v) using only v's neighbors."""
        neigh_set = {u for u in G.adj_list[v] if core_num[u] >= k}
        if len(neigh_set) < k:
            return 0
        return kCoreBaseStructuralDiversity._fast_k_core_count(
            G, neigh_set, k
        )

    # -------------------------------------------------------------------------
    # Fallback when no (k+1)-core exists in the graph
    # -------------------------------------------------------------------------
    @staticmethod
    def _fallback_k_core_method(G, k: int,
                                core_num: List[int]) -> List[int]:
        """Compute τ_k(v) directly when the graph only has a k-core."""
        n = G.vertex_num
        tau = [0] * n
        for v in range(n):
            qualified = [u for u in G.adj_list[v] if core_num[u] >= k]
            if len(qualified) >= k:
                tau[v] = kCoreBaseStructuralDiversity._compute_tau_k_local(
                    G, v, k, core_num
                )
        return tau


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
