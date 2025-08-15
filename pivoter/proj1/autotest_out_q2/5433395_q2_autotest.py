#!/usr/bin/env python3
# Auto-generated for 5433395

STUDENT_ID = "5433395"
STUDENT_NAME = "Qichen Zhang"

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
        n = G.vertex_num
        result = []

        for v in range(n):
            subG = kCoreBaseStructuralDiversity.build_neighbour_induced_subgraph(G, v)
            k_core_nodes = kCoreBaseStructuralDiversity.extract_k_core(subG.adj_list, k)
            tau_k_v = kCoreBaseStructuralDiversity.count_connected_components_union_find(subG.adj_list, k_core_nodes)
            result.append(tau_k_v)

        return result

    @staticmethod
    def build_neighbour_induced_subgraph(G, v):
        N_v = G.adj_list[v]
        neighbour_set = set(N_v)
        induced_edges = []

        # add self-loops for each vertex in N_v
        for u in N_v:
            induced_edges.append((u, u))

        for u in N_v:
            for w in G.adj_list[u]:
                if w in neighbour_set and u < w:
                    induced_edges.append((u, w))

        vertex_num = max(N_v) + 1 if N_v else v + 1
        header = (vertex_num, len(induced_edges))
        edge_list = [header] + induced_edges
        return UndirectedUnweightedGraph(edge_list)

    @staticmethod
    def extract_k_core(adj_list, k):
        from collections import deque

        n = len(adj_list)
        deg = [0] * n
        removed = [False] * n

        for u in range(n):
            deg_u = len(adj_list[u])

            deg[u] = deg_u - 2


        queue = deque([u for u in range(n) if deg[u] < k])

        while queue:
            u = queue.popleft()
            removed[u] = True
            for v in adj_list[u]:
                if v == u or removed[v]:
                    continue
                deg[v] -= 1
                if deg[v] == k - 1:
                    queue.append(v)

        return {u for u in range(n) if not removed[u]}

    @staticmethod
    def count_connected_components_union_find(adj_list, nodes):
        class UnionFind:
            def __init__(self, vertices):
                self.parent = {v: v for v in vertices}

            def find(self, x):
                if self.parent[x] != x:
                    self.parent[x] = self.find(self.parent[x])
                return self.parent[x]

            def union(self, x, y):
                px, py = self.find(x), self.find(y)
                if px != py:
                    self.parent[py] = px

            def count_components(self):
                roots = set()
                for v in self.parent:
                    roots.add(self.find(v))
                return len(roots)

        uf = UnionFind(nodes)
        for u in nodes:
            for v in adj_list[u]:
                if v == u or v not in nodes:
                    continue
                uf.union(u, v)

        return uf.count_components()


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
