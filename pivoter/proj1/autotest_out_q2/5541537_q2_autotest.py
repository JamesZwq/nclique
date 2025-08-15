#!/usr/bin/env python3
# Auto-generated for 5541537

STUDENT_ID = "5541537"
STUDENT_NAME = "Guotao Xie"

# ======= 学生代码 =======
from collections import deque
from typing import List, Set

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k: int) -> List[int]:
        n = G.vertex_num
        result = [0] * n
        if k < 0 or n == 0:
            return result

        # Global core-number decomposition (O(n + m))
        core_vals = kCoreBaseStructuralDiversity._get_core_numbers(G)

        # Collect vertices in the (k+1)-core
        high_core = {v for v, c in enumerate(core_vals) if c >= k + 1}
        if not high_core:
            # No (k+1)-core → fallback to local peeling for all vertices
            for v in range(n):
                result[v] = kCoreBaseStructuralDiversity._peel_and_count(G, v, k)
            return result

        # Find connected components in (k+1)-core subgraph
        components = kCoreBaseStructuralDiversity._get_subgraph_components(G, high_core)

        # Batch-process each component
        processed = set()
        for comp in components:
            kCoreBaseStructuralDiversity._evaluate_component(G, comp, k, result)
            processed |= comp

        # Process vertex not in those components
        for v in range(n):
            if v not in processed:
                result[v] = kCoreBaseStructuralDiversity._peel_and_count(G, v, k)

        return result

    @staticmethod
    def _get_core_numbers(G) -> List[int]:
        n = G.vertex_num
        deg = [len(G.adj_list[v]) for v in range(n)]
        max_deg = max(deg) if deg else 0

        # Bucket vertices by current degree
        buckets = [[] for _ in range(max_deg + 1)]
        for v, d in enumerate(deg):
            buckets[d].append(v)

        core = [0] * n
        # Peel vertices in non-decreasing degree order
        for d in range(max_deg + 1):
            bucket = buckets[d]
            while bucket:
                v = bucket.pop()
                core[v] = d
                for w in G.adj_list[v]:
                    if deg[w] > d:
                        old = deg[w]
                        deg[w] -= 1
                        buckets[old - 1].append(w)
        return core

    @staticmethod
    def _get_subgraph_components(G, verts: Set[int]) -> List[Set[int]]:
        visited = set()
        comps = []
        for v in verts:
            if v in visited:
                continue
            queue = deque([v])
            visited.add(v)
            comp = {v}
            while queue:
                u = queue.popleft()
                for w in G.adj_list[u]:
                    if w in verts and w not in visited:
                        visited.add(w)
                        comp.add(w)
                        queue.append(w)
            comps.append(comp)
        return comps

    @staticmethod
    def _evaluate_component(G, comp: Set[int], k: int, result: List[int]) -> None:
        # Compute t_k(v) for each v in this component
        for v in comp:
            neigh = [w for w in G.adj_list[v] if w in comp]
            if len(neigh) < k:
                result[v] = 0
                continue

            # Build degrees in the induced neighbor subgraph
            deg_map = {u: 0 for u in neigh}
            for u in neigh:
                for w in G.adj_list[u]:
                    if w in deg_map:
                        deg_map[u] += 1

            # k-core peeling within this neighbor set
            alive = set(neigh)
            queue = deque(u for u in neigh if deg_map[u] < k)
            while queue:
                u = queue.popleft()
                if u not in alive:
                    continue
                alive.remove(u)
                for w in G.adj_list[u]:
                    if w in alive:
                        deg_map[w] -= 1
                        if deg_map[w] == k - 1:
                            queue.append(w)

            # Count connected components among the survivors
            seen = set()
            count = 0
            for u in alive:
                if u in seen:
                    continue
                count += 1
                stack = [u]
                seen.add(u)
                while stack:
                    x = stack.pop()
                    for w in G.adj_list[x]:
                        if w in alive and w not in seen:
                            seen.add(w)
                            stack.append(w)
            result[v] = count

    @staticmethod
    def _peel_and_count(G, v: int, k: int) -> int:
        neigh = G.adj_list[v]
        if len(neigh) < k:
            return 0

        deg_map = {u: 0 for u in neigh}
        for u in neigh:
            for w in G.adj_list[u]:
                if w in deg_map:
                    deg_map[u] += 1

        alive = set(neigh)
        queue = deque(u for u in neigh if deg_map[u] < k)
        while queue:
            u = queue.popleft()
            if u not in alive:
                continue
            alive.remove(u)
            for w in G.adj_list[u]:
                if w in alive:
                    deg_map[w] -= 1
                    if deg_map[w] == k - 1:
                        queue.append(w)

        seen = set()
        count = 0
        for u in alive:
            if u in seen:
                continue
            count += 1
            stack = [u]
            seen.add(u)
            while stack:
                x = stack.pop()
                for w in G.adj_list[x]:
                    if w in alive and w not in seen:
                        seen.add(w)
                        stack.append(w)
        return count

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
