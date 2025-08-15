#!/usr/bin/env python3
# Auto-generated for 5254061

STUDENT_ID = "5254061"
STUDENT_NAME = "Yang Hu"

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
        n = G.vertex_num
        tau = [0] * n  # Using tau instead of sd for clarity

        # Handle edge case: empty graph
        if n == 0:
            return tau

        # For each vertex v, compute τ_k(v)
        for v in range(n):                              # O(n) iterations
            neighbors = G.adj_list[v]                   # O(1)

            # Handle isolated vertices
            if len(neighbors) == 0:
                tau[v] = 0
                continue

            # Build neighbor-induced subgraph G[N(v)]    # O(d²)
            neighbor_set = set(neighbors)
            neighbor_list = list(neighbors)
            neighbor_to_idx = {neighbor: idx for idx, neighbor in enumerate(neighbor_list)}

            # Create adjacency list for neighbor-induced subgraph
            # Use sets to avoid duplicate edges
            subgraph_adj = [set() for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list):
                for neighbor in G.adj_list[u]:
                    if neighbor in neighbor_set:
                        j = neighbor_to_idx[neighbor]
                        if i != j:  # Avoid self-loops
                            subgraph_adj[i].add(j)

            # Convert sets to lists for processing
            subgraph_adj = [list(s) for s in subgraph_adj]

            # Count k-cores in the neighbor-induced subgraph     # O(m')
            tau[v] = kCoreBaseStructuralDiversity._count_kcores(subgraph_adj, k)

        return tau

    @staticmethod
    def _count_kcores(adj_list, k):
        n = len(adj_list)
        if n == 0:
            return 0

        # Special case: k <= 0, every connected component is a k-core
        if k <= 0:
            return kCoreBaseStructuralDiversity._count_all_connected_components(adj_list)

        # Step 1: Find all vertices that can be in some k-core using core decomposition
        core_numbers = kCoreBaseStructuralDiversity._compute_core_numbers(adj_list)                   # O(m')

        # Step 2: Get vertices with core number >= k
        kcore_candidates = {i for i in range(n) if core_numbers[i] >= k}

        if not kcore_candidates:
            return 0

        # Step 3: Find connected components among candidates and verify they form valid k-cores
        return kCoreBaseStructuralDiversity._count_valid_kcores(adj_list, kcore_candidates, k)       # O(d + m')

    @staticmethod
    def _compute_core_numbers(adj_list):

        n = len(adj_list)
        degrees = [len(adj_list[i]) for i in range(n)]
        core_numbers = [0] * n
        processed = [False] * n

        # Create buckets for vertices by degree
        max_degree = max(degrees) if degrees else 0
        buckets = [[] for _ in range(max_degree + 1)]    # O(d)

        for i in range(n):
            buckets[degrees[i]].append(i)

        # Process vertices in order of increasing degree
        for deg in range(max_degree + 1):                # O(d + m')
            while buckets[deg]:
                v = buckets[deg].pop()
                if processed[v]:
                    continue

                core_numbers[v] = degrees[v]
                processed[v] = True

                # Update degrees of neighbors
                for u in adj_list[v]:
                    if not processed[u] and degrees[u] > degrees[v]:
                        # Remove u from its current bucket
                        old_deg = degrees[u]
                        degrees[u] -= 1
                        new_deg = degrees[u]

                        if new_deg >= deg:
                            buckets[new_deg].append(u)

        return core_numbers

    @staticmethod
    def _count_valid_kcores(adj_list, candidates, k):
        if not candidates:
            return 0

        visited = set()
        kcore_count = 0

        for start in candidates:
            if start not in visited:
                # BFS to find connected component
                component = set()
                queue = deque([start])
                visited.add(start)
                component.add(start)

                while queue:                     # O(d + m')
                    u = queue.popleft()
                    for v in adj_list[u]:
                        if v in candidates and v not in visited:
                            visited.add(v)
                            component.add(v)
                            queue.append(v)

                # Check if this component forms a valid k-core
                # (all vertices in component should have degree >= k within the component)
                if kCoreBaseStructuralDiversity._is_valid_kcore(adj_list, component, k):      # O(d²) worst case
                    kcore_count += 1

        return kcore_count

    @staticmethod
    def _is_valid_kcore(adj_list, component, k):
        for v in component:
            degree_in_component = sum(1 for u in adj_list[v] if u in component)
            if degree_in_component < k:
                return False
        return True

    @staticmethod
    def _count_all_connected_components(adj_list):

        n = len(adj_list)
        if n == 0:
            return 0

        visited = [False] * n
        components = 0

        for start in range(n):
            if not visited[start]:
                # BFS to find connected component
                queue = deque([start])
                visited[start] = True

                while queue:
                    u = queue.popleft()
                    for v in adj_list[u]:
                        if not visited[v]:
                            visited[v] = True
                            queue.append(v)

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
