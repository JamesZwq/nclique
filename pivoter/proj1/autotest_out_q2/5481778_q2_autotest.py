#!/usr/bin/env python3
# Auto-generated for 5481778

STUDENT_ID = "5481778"
STUDENT_NAME = "Taylor Hoffman"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # u refers to some vertex's neighbour
    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n
        for v in range(n):
          sd[v] = kCoreBaseStructuralDiversity.k_core_base_components(G, v, k, n)
        return sd

    @staticmethod
    def k_core_base_components(G, curr_vertex, k, n):

        neighbours = set(G.adj_list[curr_vertex])
        # No neighbours
        if not neighbours:
            return 0

        # Create subgraph of neighbours using adj_list and degrees
        subgraph_nbrs = [[] for _ in range(n)]
        subgraph_nbrs_degrees = [0] * n

        for nbr in neighbours:
            for u in G.adj_list[nbr]:
                # Create edge (nbr, u) only if u is a neighbour and u > nbr to avoid duplicates
                if u in neighbours and u > nbr:
                    subgraph_nbrs[nbr].append(u)
                    subgraph_nbrs_degrees[nbr] += 1

                    subgraph_nbrs[u].append(nbr)
                    subgraph_nbrs_degrees[u] += 1

        # Create initial queue for all vertices with degree < k
        dq = deque()
        visited = set()
        for nbr in neighbours:
          if subgraph_nbrs_degrees[nbr] < k:
            dq.append(nbr)
            visited.add(nbr)

        # Peel all vertices with degree < k by 1
        while dq:
          curr_nbr = dq.popleft()
          for u in subgraph_nbrs[curr_nbr]:
            if u not in visited:
              subgraph_nbrs_degrees[u] -= 1
              if subgraph_nbrs_degrees[u] < k:
                dq.append(u)
                visited.add(u)

        # Count all nbr connected components with at least a k-core
        visited = set()
        k_core_base_components_count = 0
        for nbr in neighbours:
            if nbr not in visited and subgraph_nbrs_degrees[nbr] >= k:
                nbr_component = set()
                kCoreBaseStructuralDiversity.explore_k_core_subgraph(G, nbr, k, subgraph_nbrs, subgraph_nbrs_degrees, nbr_component, visited)
                k_core_base_components_count += 1

        return k_core_base_components_count

    @staticmethod
    def explore_k_core_subgraph(G, nbr, k, subgraph_nbrs, subgraph_nbrs_degrees, nbr_component, visited):
        # Explore using dfs to fina a single k-core component
        nbr_component.add(nbr)
        visited.add(nbr)
        for u in subgraph_nbrs[nbr]:
            if u not in visited and subgraph_nbrs_degrees[u] >= k:
                kCoreBaseStructuralDiversity.explore_k_core_subgraph(G, u, k, subgraph_nbrs, subgraph_nbrs_degrees, nbr_component, visited)

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
