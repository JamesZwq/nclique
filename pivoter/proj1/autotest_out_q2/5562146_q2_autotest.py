#!/usr/bin/env python3
# Auto-generated for 5562146

STUDENT_ID = "5562146"
STUDENT_NAME = "Zhengtao Ying"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core-based structural diversity for each vertex.
        Parameters:
        -----------
        G : UndirectedUnweightedGraph
            The input undirected, unweighted graph.
        k : int
            The minimum degree for core definition.
        Returns:
        --------
        List[int]
            A list where the i-th element is τ_k(i), the number of k-cores
            in the neighbor-induced subgraph of vertex i.
        """

        def neighbor_subgraph(v):
            neighbors = G.adj_list[v]
            edge_map = defaultdict(set)
            neighbor_set = set(neighbors)
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbor_set and u < w:
                        edge_map[u].add(w)
                        edge_map[w].add(u)
            return neighbors, edge_map

        def count_k_cores(nodes, edges):
            if not nodes:
                return 0

            degree = {node: len(edges[node]) for node in nodes}
            alive = set(nodes)
            queue = deque([v for v in nodes if degree[v] < k])

            while queue:
                u = queue.popleft()
                if u not in alive:
                    continue
                alive.remove(u)
                for nei in edges[u]:
                    if nei in alive:
                        degree[nei] -= 1
                        if degree[nei] < k:
                            queue.append(nei)

            return connected_components(alive, edges) if alive else 0

        def connected_components(nodes, edges):
            seen = set()
            components = 0

            for v in nodes:
                if v not in seen:
                    components += 1
                    queue = deque([v])
                    seen.add(v)
                    while queue:
                        cur = queue.popleft()
                        for nbr in edges[cur]:
                            if nbr in nodes and nbr not in seen:
                                seen.add(nbr)
                                queue.append(nbr)
            return components

        result = [0] * G.vertex_num
        for v in range(G.vertex_num):
            neighbors, edges = neighbor_subgraph(v)
            result[v] = count_k_cores(neighbors, edges)
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
