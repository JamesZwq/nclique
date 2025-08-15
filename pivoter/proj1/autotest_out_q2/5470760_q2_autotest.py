#!/usr/bin/env python3
# Auto-generated for 5470760

STUDENT_ID = "5470760"
STUDENT_NAME = "Yunpu Du"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################

from collections import deque, defaultdict

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n
        if k < 0 or n == 0:
            return sd

        core_numbers = kCoreBaseStructuralDiversity._compute_core_numbers(G)
        k_plus_1_vertices = set(v for v in range(n) if core_numbers[v] >= k + 1)

        if not k_plus_1_vertices:
            return kCoreBaseStructuralDiversity._fallback_k_core_method(G, k, core_numbers)

        k_plus_1_components = kCoreBaseStructuralDiversity._generic_traverse_components(
            G, k_plus_1_vertices, collect=True)

        for component in k_plus_1_components:
            kCoreBaseStructuralDiversity._process_component_batch(G, component, k, sd)

        processed = {u for comp in k_plus_1_components for v in comp for u in G.adj_list[v]}

        for v in range(n):
            if v not in processed or sd[v] == 0:
                sd[v] += kCoreBaseStructuralDiversity._compute_tau_k_local(G, v, k, core_numbers)

        return sd

    @staticmethod
    def _compute_core_numbers(G):
        n = G.vertex_num
        degree = [len(G.adj_list[v]) for v in range(n)]
        core_number = [0] * n
        max_degree = max(degree, default=0)

        count = [0] * (max_degree + 1)
        for d in degree:
            count[d] += 1

        start = [0] * (max_degree + 1)
        for i in range(1, max_degree + 1):
            start[i] = start[i - 1] + count[i - 1]

        sorted_vertices = [0] * n
        pos = start[:]
        for v in range(n):
            sorted_vertices[pos[degree[v]]] = v
            pos[degree[v]] += 1

        for v in sorted_vertices:
            core_number[v] = degree[v]
            for u in G.adj_list[v]:
                if degree[u] > degree[v]:
                    degree[u] -= 1

        return core_number

    @staticmethod
    def _process_component_batch(G, component, k, sd):
        comp_set = set(component)
        all_neighbors = {u for v in component for u in G.adj_list[v]}

        for v in all_neighbors:
            neighbors = set(G.adj_list[v])
            component_neighbors = neighbors & comp_set
            if len(component_neighbors) >= k:
                k_core_count = kCoreBaseStructuralDiversity._fast_k_core_count(
                    G, component_neighbors, k, neighbors)
                sd[v] += k_core_count

    @staticmethod
    def _fast_k_core_count(G, candidates, k, valid_neighbors):
        if len(candidates) < k:
            return 0

        for v in candidates:
            deg = sum(1 for u in G.adj_list[v] if u in candidates and u in valid_neighbors)
            if deg < k:
                return kCoreBaseStructuralDiversity._exact_k_core_count(G, candidates, k, valid_neighbors)

        return kCoreBaseStructuralDiversity._generic_traverse_components(
            G, candidates, valid_neighbors, collect=False)

    @staticmethod
    def _exact_k_core_count(G, vertices, k, valid_neighbors):
        k_core_vertices = kCoreBaseStructuralDiversity._peeling_k_core(vertices, G.adj_list, k, valid_neighbors)
        if not k_core_vertices:
            return 0

        return kCoreBaseStructuralDiversity._generic_traverse_components(
            G, k_core_vertices, valid_neighbors, collect=False)

    @staticmethod
    def _compute_tau_k_local(G, v, k, core_numbers):
        neighbors = [u for u in G.adj_list[v] if core_numbers[u] >= k]
        if len(neighbors) < k:
            return 0

        neighbor_set = set(neighbors)
        k_core_vertices = kCoreBaseStructuralDiversity._peeling_k_core(neighbor_set, G.adj_list, k)
        if not k_core_vertices:
            return 0

        return kCoreBaseStructuralDiversity._generic_traverse_components(
            G, k_core_vertices, collect=False)

    @staticmethod
    def _fallback_k_core_method(G, k, core_numbers):
        sd = [0] * G.vertex_num
        for v in range(G.vertex_num):
            neighbors = [u for u in G.adj_list[v] if core_numbers[u] >= k]
            if len(neighbors) >= k:
                sd[v] = kCoreBaseStructuralDiversity._compute_tau_k_local(G, v, k, core_numbers)
        return sd

    @staticmethod
    def _peeling_k_core(vertices, adj_list, k, valid_neighbors=None):
        degree = {
            v: sum(1 for u in adj_list[v] if u in vertices and (valid_neighbors is None or u in valid_neighbors))
            for v in vertices
        }
        queue = deque(v for v in vertices if degree[v] < k)
        removed = set(queue)

        while queue:
            v = queue.popleft()
            for u in adj_list[v]:
                if u in vertices and u not in removed and (valid_neighbors is None or u in valid_neighbors):
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        removed.add(u)

        return vertices - removed

    @staticmethod
    def _generic_traverse_components(G, vertices, valid_neighbors=None, collect=False):
        visited = set()
        result = [] if collect else 0
        for start in vertices:
            if start not in visited:
                container = deque([start])
                visited.add(start)
                component = [start] if collect else None

                while container:
                    u = container.popleft()
                    for v in G.adj_list[u]:
                        if v in vertices and v not in visited and (valid_neighbors is None or v in valid_neighbors):
                            visited.add(v)
                            container.append(v)
                            if collect:
                                component.append(v)

                if collect:
                    result.append(component)
                else:
                    result += 1

        return result


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
