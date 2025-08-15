#!/usr/bin/env python3
# Auto-generated for 5468081

STUDENT_ID = "5468081"
STUDENT_NAME = "Tinghao Li"

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
        num_vertices = G.vertex_num
        results = [0] * num_vertices

        for node in range(num_vertices):
            if not G.adj_list[node]:
                continue

            subgraph = kCoreBaseStructuralDiversity._neighbor_induced_subgraph(G, node)
            core_count = kCoreBaseStructuralDiversity._count_kcore_components(subgraph, k)
            results[node] = core_count

        return results

    @staticmethod
    def _neighbor_induced_subgraph(G, center):
        """Build the subgraph induced by the neighbors of a vertex."""
        neighbors = set(G.adj_list[center])
        induced = defaultdict(set)

        for u in neighbors:
            for v in G.adj_list[u]:
                if v in neighbors:
                    induced[u].add(v)

        return induced

    @staticmethod
    def _count_kcore_components(adj_list, k):
        """Count the number of connected components in the k-core of the given subgraph."""
        if not adj_list:
            return 0

        degrees = {node: len(neighs) for node, neighs in adj_list.items()}
        to_remove = deque([v for v, d in degrees.items() if d < k])
        removed = set(to_remove)

        # Iteratively remove low-degree nodes
        while to_remove:
            curr = to_remove.popleft()
            for neighbor in adj_list[curr]:
                if neighbor not in removed:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        removed.add(neighbor)
                        to_remove.append(neighbor)

        remaining = set(adj_list.keys()) - removed
        if not remaining:
            return 0

        # Count connected components in remaining subgraph
        visited = set()
        components = 0

        for node in remaining:
            if node in visited:
                continue
            components += 1
            stack = [node]
            visited.add(node)

            while stack:
                u = stack.pop()
                for v in adj_list[u]:
                    if v not in removed and v not in visited:
                        visited.add(v)
                        stack.append(v)

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
