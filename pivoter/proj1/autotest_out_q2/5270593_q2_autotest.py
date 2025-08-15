#!/usr/bin/env python3
# Auto-generated for 5270593

STUDENT_ID = "5270593"
STUDENT_NAME = "Mingyuan Cui"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque

################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def k_core(G: dict[int, list[int]], k: int) -> dict[int, list[int]]:
        deg = {v: len(G[v]) for v in G}  # O(|V|)
        queue = deque([v for v in G.keys() if deg[v] < k])  # O(|V|)
        removed = {v: False for v in G.keys()}  # O(|V|)
        for v in queue:
            removed[v] = True  # O(|V|)
        while queue:  # At most |V| iterations total
            v = queue.popleft()  # O(1)
            for u in G[v]:  # Each edge is traversed at most once, total O(|E|)
                if not removed[u]:
                    deg[u] -= 1
                    if deg[u] < k:
                        queue.append(u)
                        removed[u] = True
        k_core_graph = {v: [] for v in G.keys() if not removed[v]}  # O(|V|)
        for v in k_core_graph.keys():  # O(|V|)
            for u in G[v]:  # O(|E|)
                if not removed[u]:
                    k_core_graph[v].append(u)
        return k_core_graph

    @staticmethod
    def subgraph(G, vertices: list[int]) -> dict[int, list[int]]:
        vertex_set = set(vertices)  # O(|V'|)
        subgraph_adj = {v: [] for v in vertices}  # O(|V'|)
        for u in vertex_set:  # O(|V'|)
            for v in G.adj_list[u]:  # Total O(|E'|)
                if v in vertex_set:  # O(1) - set lookup
                    subgraph_adj[u].append(v)
        return subgraph_adj

    @staticmethod
    def parallel_connected_components(G: dict[int, list[int]]) -> list[list[int]]:
        import concurrent.futures
        from collections import deque
        unique_components = set()
        components = []
        def find_component_bfs(start_v):
            local_visited = set([start_v])
            component = [start_v]
            queue = deque([start_v])
            while queue:
                v = queue.popleft()
                for u in G[v]:
                    if u not in local_visited:
                        local_visited.add(u)
                        component.append(u)
                        queue.append(u)
            return component, local_visited
        visited = set()
        unvisited_nodes = list(G.keys())
        with concurrent.futures.ThreadPoolExecutor() as executor:
            while unvisited_nodes:
                candidates = []
                for v in list(unvisited_nodes):
                    if v not in visited:
                        candidates.append(v)
                        unvisited_nodes.remove(v)
                        if len(candidates) >= min(8, len(G)):
                            break
                if not candidates:
                    break
                future_to_node = {executor.submit(find_component_bfs, v): v for v in candidates}
                for future in concurrent.futures.as_completed(future_to_node):
                    component, local_visited = future.result()
                    component_set = frozenset(component)
                    if component_set not in unique_components:
                        unique_components.add(component_set)
                        components.append(component)
                    visited.update(local_visited)
                    unvisited_nodes = [v for v in unvisited_nodes if v not in visited]
        return components

    @staticmethod
    def connected_components(G: dict[int, list[int]]) -> list[list[int]]:
        visited: set[int] = set()
        components: list[list[int]] = []

        def dfs(v: int, component: list[int]):
            visited.add(v)
            component.append(v)
            for u in G[v]:
                if u not in visited:
                    dfs(u, component)

        for v in G.keys():  # O(|V|)
            if v not in visited:
                component = []
                dfs(v, component)  # Total DFS complexity is O(|V|+|E|)
                components.append(component)
        return components # O(|V|+|E|)

    @staticmethod
    def process(G, k: int) -> list[int]:
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int

        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        n = G.vertex_num
        sd = [0] * n  # O(|V|)
        for v in range(0, n):  # O(|V|)
            nbrs = G.adj_list[v]
            # Extract subgraph - O(|V|+|E|)
            subgraph = kCoreBaseStructuralDiversity.subgraph(G, nbrs)
            # Compute k-core - O(|V|+|E|)
            k_core = kCoreBaseStructuralDiversity.k_core(subgraph, k)
            # Find connected components - O(|V|+|E|)
            if k > 10:
                components = kCoreBaseStructuralDiversity.parallel_connected_components(k_core)
            else:
                components = kCoreBaseStructuralDiversity.connected_components(k_core)
            count = len(components)
            sd[v] = count
        return sd  # O(|V| * (|V| + |E|)) = O(|V|^2 + |V| * |E|), worst case O(|V|^2)

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
