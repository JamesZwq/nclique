#!/usr/bin/env python3
# Auto-generated for 5484763

STUDENT_ID = "5484763"
STUDENT_NAME = "Jinrui Cheng"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    _neighbor_cache = {}
    _current_graph_id = None
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
        # Clear cache when graph changes
        graph_id = id(G)
        if kCoreBaseStructuralDiversity._current_graph_id != graph_id:
            kCoreBaseStructuralDiversity._neighbor_cache.clear()
            kCoreBaseStructuralDiversity._current_graph_id = graph_id

        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]

            # Optimization: directly check empty neighbors
            if not neighbors:
                continue  # sd[v] is already 0

            # Cache key: sorted tuple of neighbors
            cache_key = tuple(sorted(neighbors))

            if cache_key not in kCoreBaseStructuralDiversity._neighbor_cache:
                neighbor_graph = kCoreBaseStructuralDiversity._build_neighbor_induced_subgraph(G, neighbors)
                kCoreBaseStructuralDiversity._neighbor_cache[cache_key] = neighbor_graph
            else:
                neighbor_graph = kCoreBaseStructuralDiversity._neighbor_cache[cache_key]

            # Count k-cores
            sd[v] = kCoreBaseStructuralDiversity._count_k_cores(neighbor_graph, k)

        return sd

    @staticmethod
    def _build_neighbor_induced_subgraph(G, neighbors):
        if not neighbors:
            return {}

        # Use defaultdict to automatically create empty lists, avoiding manual initialization
        nbr_graph = defaultdict(list)
        neighbor_set = set(neighbors)

        # Build adjacency list of neighbor-induced subgraph (keeping original vertex IDs)
        for u in neighbors:
            for w in G.adj_list[u]:
                # If w is also a neighbor and not u itself
                if w in neighbor_set and w != u:
                    nbr_graph[u].append(w)

        # Convert to regular dict (ensure all neighbors have keys)
        result = {}
        for u in neighbors:
            result[u] = nbr_graph[u]

        return result

    @staticmethod
    def _count_k_cores(nbr_graph, k):
        """
        Count the number of k-cores (connected components)

        Parameters
        ----------
        nbr_graph : Dict[int, List[int]] - adjacency list of the graph
        k : int - k value

        Returns
        -------
        int - number of k-core connected components
        """
        if not nbr_graph:
            return 0

        vertices = set(nbr_graph.keys())

        # Special case for k=0: every connected component is a 0-core
        if k == 0:
            return kCoreBaseStructuralDiversity._count_connected_components(nbr_graph)

        # Calculate initial degrees
        degrees = {u: len(nbr_graph[u]) for u in vertices}

        # k-core decomposition: iteratively remove vertices with degree < k
        removed = set()
        queue = deque()

        # Initialize: add all vertices with degree < k to removal queue
        for u in vertices:
            if degrees[u] < k:
                queue.append(u)
                removed.add(u)

        # Perform k-core decomposition
        while queue:
            u = queue.popleft()
            # Update degrees of all neighbors of u
            for v in nbr_graph[u]:
                if v not in removed:
                    degrees[v] -= 1
                    # If neighbor's degree drops below k, also remove it
                    if degrees[v] < k:
                        queue.append(v)
                        removed.add(v)

        # Optimization: use set difference to get remaining vertices
        remaining = vertices - removed
        if not remaining:
            return 0

        # Build subgraph of remaining vertices (using basic loops to avoid complex dict comprehension)
        remaining_graph = {}
        for u in remaining:
            remaining_graph[u] = []
            for v in nbr_graph[u]:
                if v in remaining:
                    remaining_graph[u].append(v)

        return kCoreBaseStructuralDiversity._count_connected_components(remaining_graph)

    @staticmethod
    def _count_connected_components(graph):
        if not graph:
            return 0

        unvisited = set(graph.keys())
        component_count = 0

        while unvisited:
            # Start BFS from any unvisited vertex
            start = unvisited.pop()
            queue = deque([start])

            while queue:
                u = queue.popleft()
                for v in graph[u]:
                    if v in unvisited:
                        unvisited.remove(v)
                        queue.append(v)

            component_count += 1

        return component_count


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
