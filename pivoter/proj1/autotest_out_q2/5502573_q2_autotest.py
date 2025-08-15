#!/usr/bin/env python3
# Auto-generated for 5502573

STUDENT_ID = "5502573"
STUDENT_NAME = "Wenao Yao"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

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
        # Step 1: compute core numbers for all vertices
        core_numbers = kCoreBaseStructuralDiversity._compute_core_numbers(G)

        # Step 2: collect all vertices with core number >= k + 1
        k_plus_1_vertices = set(v for v in range(n) if core_numbers[v] >= k + 1)

        # Step 3: find connected components in (k+1)-core subgraph
        k_plus_1_components = kCoreBaseStructuralDiversity._find_components(G, k_plus_1_vertices)

        # Step 4: batch process (k+1)-core components
        for component in k_plus_1_components:
            kCoreBaseStructuralDiversity._process_component_batch(G, component, k, sd)

        # Step 5: fallback for remaining vertices
        processed = set()
        for comp in k_plus_1_components:
            for v in comp:
                for u in G.adj_list[v]:
                    processed.add(u)

        for v in range(n):
            if v not in processed or sd[v] == 0:
                sd[v] += kCoreBaseStructuralDiversity._compute_tau_k_local(G, v, k, core_numbers)

        return sd

    @staticmethod
    def _compute_core_numbers(G):
        n = G.vertex_num
        if n == 0:
            return []

        degree = [len(G.adj_list[v]) for v in range(n)]
        core_number = [0] * n

        max_degree = max(degree) if degree else 0
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

        for i in range(n):
            v = sorted_vertices[i]
            core_number[v] = degree[v]
            for u in G.adj_list[v]:
                if degree[u] > degree[v]:
                    degree[u] -= 1

        return core_number

    @staticmethod
    def _find_components(G, vertices):
        visited = set()
        components = []
        for v in vertices:
            if v not in visited:
                component = []
                queue = deque([v])
                visited.add(v)
                while queue:
                    u = queue.popleft()
                    component.append(u)
                    for w in G.adj_list[u]:
                        if w in vertices and w not in visited:
                            queue.append(w)
                            visited.add(w)
                components.append(component)
        return components

    @staticmethod
    def _process_component_batch(G, component, k, sd):
        component_set = set(component)
        all_neighbors = set()
        for v in component:
            all_neighbors.update(G.adj_list[v])

        for v in all_neighbors:
            neighbors = set(G.adj_list[v])
            component_neighbors = neighbors & component_set
            if len(component_neighbors) >= k:
                k_core_count = kCoreBaseStructuralDiversity._fast_k_core_count(
                    G, component_neighbors, k, neighbors
                )
                sd[v] += k_core_count

    @staticmethod
    def _fast_k_core_count(G, candidates, k, valid_neighbors):
        if len(candidates) < k:
            return 0

        subgraph_degree = {}
        for v in candidates:
            degree = sum(1 for u in G.adj_list[v] if u in candidates and u in valid_neighbors)
            if degree < k:
                return kCoreBaseStructuralDiversity._exact_k_core_count(
                    G, candidates, k, valid_neighbors
                )
            subgraph_degree[v] = degree

        return kCoreBaseStructuralDiversity._count_components(G, candidates, valid_neighbors)

    @staticmethod
    def _exact_k_core_count(G, vertices, k, valid_neighbors):
        vertices = set(vertices)
        degree = {}
        for v in vertices:
            degree[v] = sum(1 for u in G.adj_list[v] if u in vertices and u in valid_neighbors)

        queue = deque(v for v in vertices if degree[v] < k)
        removed = set(queue)

        while queue:
            v = queue.popleft()
            for u in G.adj_list[v]:
                if u in vertices and u not in removed and u in valid_neighbors:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        removed.add(u)

        k_core_vertices = vertices - removed
        if not k_core_vertices:
            return 0

        return kCoreBaseStructuralDiversity._count_components(G, k_core_vertices, valid_neighbors)

    @staticmethod
    def _count_components(G, vertices, valid_neighbors):
        visited = set()
        components = 0
        for v in vertices:
            if v not in visited:
                queue = deque([v])
                visited.add(v)
                while queue:
                    u = queue.popleft()
                    for w in G.adj_list[u]:
                        if w in vertices and w in valid_neighbors and w not in visited:
                            queue.append(w)
                            visited.add(w)
                components += 1
        return components

    @staticmethod
    def _compute_tau_k_local(G, v, k, core_numbers):
        neighbors = set(G.adj_list[v])
        if not neighbors:
            return 0

        nbr_graph = {u: [] for u in neighbors}
        for u in neighbors:
            for w in G.adj_list[u]:
                if w in neighbors:
                    nbr_graph[u].append(w)

        k_core_vertices = kCoreBaseStructuralDiversity._compute_k_core(nbr_graph, k)
        if not k_core_vertices:
            return 0

        return len(kCoreBaseStructuralDiversity._find_connected_components(nbr_graph, k_core_vertices))

    @staticmethod
    def _compute_k_core(graph, k):
        degree  = {v: len(adj) for v, adj in graph.items()}
        queue   = deque([v for v, deg in degree.items() if deg < k])
        removed = set(queue)

        while queue:
            v = queue.popleft()
            for u in graph[v]:
                if u not in removed:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        removed.add(u)

        return set(graph.keys()) - removed

    @staticmethod
    def _find_connected_components(graph: dict, vertices: set):
        visited = set()
        components = []
        for v in vertices:
            if v not in visited:
                comp = set()
                queue = deque([v])
                visited.add(v)
                while queue:
                    u = queue.popleft()
                    comp.add(u)
                    for w in graph[u]:
                        if w in vertices and w not in visited:
                            visited.add(w)
                            queue.append(w)
                components.append(comp)
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
