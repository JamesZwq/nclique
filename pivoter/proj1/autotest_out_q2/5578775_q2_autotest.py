#!/usr/bin/env python3
# Auto-generated for 5578775

STUDENT_ID = "5578775"
STUDENT_NAME = "Marcel Juventin"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
from typing import DefaultDict, Set
from dataclasses import dataclass
################################################################################

class kCoreBaseStructuralDiversity(object):

    class kCoreDecomposition():
        @dataclass
        class State:
            """Corresponds to data structures outlined in Khaouid et al. 2015"""
            d: list[int]
            b: list[int]
            D: list[int]
            p: list[int]

        def __init__(self, G, k: int, vertex: int):
            self.G = G
            self.k = k
            self.n = G.vertex_num
            self.vertex = vertex
            self.label_to_index: dict[int] = {}
            self.index_to_label: dict[int] = {}
            self.nbr_adj_list: defaultdict[int, list[int]] = defaultdict(list)

            self.state = self.initialize_structures()


        def set_bidirectional_map(self, nbr_map: list[int]) -> None:
            """Helper method to map vertex Ids to zero based indexing"""
            self.label_to_index = {v: i for i, v in enumerate(nbr_map)}
            self.index_to_label = {i: v for i, v in enumerate(nbr_map)}


        def initialize_structures(self) -> None:
            """Initialization step"""

            n = self.n
            adj_list = self.G.adj_list

            neighbour_list = adj_list[self.vertex] if adj_list else []
            n_neighbours = len(neighbour_list)

            # Handle edge cases, empty adjacency list or no edges.
            if n_neighbours == 0 or not adj_list:
                return self.State([], [], [], [])

            D_ptr = 0
            in_neighbourhood = [False] *  n
            queue = deque()
            nbr_map = [None] * n_neighbours

            p = [-1] *  n_neighbours
            b = [-1] * n_neighbours

            D = [None] * n_neighbours
            d = [0] *  n_neighbours


            # Obtains nbr_v and populates structure to map vertex id to index
            for v in neighbour_list:#self.G.adj_list[self.vertex]:
                in_neighbourhood[v] = True
                nbr_map[D_ptr] = v
                D[D_ptr] = D_ptr
                p[D_ptr] = D_ptr
                D_ptr += 1

            # Create lookup dicts bidirectionally id -> index and vice versa.
            self.set_bidirectional_map(nbr_map)

            # Add edges and increment degree to subgraph iff u and v are neighbours of vertex.
            for u in neighbour_list:
                for v in self.G.adj_list[u]:
                    if in_neighbourhood[v] and u < v:
                        self.nbr_adj_list[self.label_to_index[u]].append(self.label_to_index[v])
                        self.nbr_adj_list[self.label_to_index[v]].append(self.label_to_index[u])
                        d[self.label_to_index[u]] += 1
                        d[self.label_to_index[v]] += 1

            # Obtain maximum degree, set up bin offsets.
            d_max = max(d)
            b = [0] * (d_max + 1)
            for x in d:
                b[x] += 1

            # Use bin offsets to intialize bin indices.
            idx = 0
            for i in range(d_max + 1):
                offset = b[i]
                b[i] = idx
                idx += offset

            # Create auxiliary bin use it to arrange vertices and positions in correct order.
            aux_b = b[:]
            for i in range(n_neighbours):
                di = d[i]
                idx = aux_b[di]
                D[idx] = i
                p[i] = idx
                aux_b[di] += 1

            return self.State(d, b, D, p)


        def perform_decomposition(self) -> list[int]:
            """Implementation of algorithm introduced by Khaouid et al. 2015"""
            D = self.state.D
            d = self.state.d
            p = self.state.p
            b = self.state.b
            adj_list = self.nbr_adj_list
            n_vertex_neighbours = len(adj_list)

            for i in range(n_vertex_neighbours):
                v = D[i]
                for u in adj_list[v]:
                    if d[u] > d[v]:
                        du = d[u]
                        pu = p[u]
                        pw = b[du]
                        w = D[pw]
                        if u != w:
                            D[pu] = w
                            D[pw] = u
                            p[u] = pw
                            p[w] = pu
                        b[du] += 1
                        d[u] -= 1
            return d


        def count_components(self, d: list[int]):
            """A basic BFS to count the number of k-cores."""
            adj_list = self.nbr_adj_list
            nbr_cardinality = len(adj_list)
            included = {i for i, v in enumerate(d) if v >= self.k}
            visited = [False] * nbr_cardinality
            queue = deque()

            count = 0
            for i in included:
                if not visited[i]:
                    visited[i] = True
                    queue.append(i)
                    while queue:
                        node = queue.popleft()
                        for v in adj_list[node]:
                            if not visited[v] and v in included:
                                visited[v] = True
                                queue.append(v)
                    count += 1
            return count


    def __init__(self):
        pass

    @staticmethod
    def process(G, k: int):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """

        counts = []
        for v in range(G.vertex_num):
            obj = kCoreBaseStructuralDiversity().kCoreDecomposition(G, k, v)
            d = obj.perform_decomposition()
            counts.append(obj.count_components(d))
        return counts



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
