#!/usr/bin/env python3
# Auto-generated for 5466081

STUDENT_ID = "5466081"
STUDENT_NAME = "Silin Zhang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
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
        n = G.vertex_num
        result = [0] * n
        # Calculate τ_k(v) for each vertex
        for v in range(n):
            neighbors = G.adj_list[v]
            if len(neighbors) == 0:
                result[v] = 0
                continue

            # Build neighbor-induced subgraph
            subgraph = kCoreBaseStructuralDiversity.neighbor_subgraph(G, neighbors)
            # Find k-core vertices
            k_core = kCoreBaseStructuralDiversity.find_k_core(subgraph, k)
            # Count connected components
            result[v] = kCoreBaseStructuralDiversity.count_components(subgraph, k_core)
        return result

    @staticmethod
    def neighbor_subgraph(G, neighbors):
        # Build subgraph induced by neighbors
        neighbor_set = set(neighbors)
        subgraph = defaultdict(list)
        # Add edges between neighbors
        for u in neighbor_set:
            for v in G.adj_list[u]:
                if v in neighbor_set:
                    subgraph[u].append(v)
        return subgraph

    @staticmethod
    def find_k_core(subgraph, k):
        # Find vertices that belong to k-core
        if not subgraph:
            return set()
        # Calculate initial degrees
        degrees = defaultdict(int)
        for u in subgraph:
            degrees[u] = len(subgraph[u])

        # Track active vertices
        active = set(subgraph.keys())
        changed = True
        while changed:
            changed = False
            vertices_to_remove = []
            # Find vertices with degree < k
            for u in active:
                if degrees[u] < k:
                    vertices_to_remove.append(u)
            # Remove vertices and update degrees
            for u in vertices_to_remove:
                active.remove(u)
                changed = True
                for v in subgraph[u]:
                    if v in active:
                        degrees[v] -= 1
        return active

    @staticmethod
    def count_components(subgraph, k_core):
        # Count connected components in k-core
        if not k_core:
            return 0
        visited = set()
        component_count = 0
        for start_vertex in k_core:
            if start_vertex in visited:
                continue

            # Start new component
            component_count += 1
            queue = deque()
            queue.append(start_vertex)
            visited.add(start_vertex)
            # Find all vertices in this component using BFS
            while queue:
                current = queue.popleft()
                for neighbor in subgraph[current]:
                    if neighbor in k_core and neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
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
