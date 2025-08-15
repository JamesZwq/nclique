#!/usr/bin/env python3
# Auto-generated for 5587914

STUDENT_ID = "5587914"
STUDENT_NAME = "Muqun Nan"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def subgraph(G, vertex):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        vertex : int
        Returns
        -------
        graph : Dict[int, Set[int]]
        numVertex : int
        """
        neighbors = set(G.adj_list[vertex])
        graph = {}

        for v in neighbors:
            graph[v] = set()
            for nbr in G.adj_list[v]:
                if nbr in neighbors:
                    graph[v].add(nbr)

        return graph, len(graph)

    @staticmethod
    def kCore(nVerts, graph, k):
        """
        Parameters
        ----------
        nVerts : int
        graph : Dict[int, Set[int]]
        k : int
        Returns
        -------
        Dict[int, Set[int]]
        int
        Set[int]
        """
        degree = {v: len(nbrs) for v, nbrs in graph.items()}
        deleted = set()
        queue = deque(v for v in graph if degree[v] < k)

        while queue:
            v = queue.popleft()
            if v in deleted:
                continue
            deleted.add(v)
            for nbr in graph[v]:
                if nbr not in deleted:
                    degree[nbr] -= 1
                    if degree[nbr] < k:
                        queue.append(nbr)
            nVerts -= 1

        return graph, nVerts, deleted

    @staticmethod
    def nKCores(graph, nVerts, deleted):
        """
        Parameters
        ----------
        graph: Dict[int, Set[int]]
        nVerts: int
        Returns
        -------
        int
        """
        parent = {v: -1 for v in graph if v not in deleted}
        size = {v: 1 for v in graph if v not in deleted}

        def find(v):
            if parent[v] == -1:
                return v
            parent[v] = find(parent[v])
            return parent[v]

        def union(x, y):
            xRoot = find(x)
            yRoot = find(y)
            if xRoot != yRoot:
                if size[xRoot] < size[yRoot]:
                    parent[xRoot] = yRoot
                    size[yRoot] += size[xRoot]
                else:
                    parent[yRoot] = xRoot
                    size[xRoot] += size[yRoot]

        for v in graph:
            if v in deleted:
                continue
            for nbr in graph[v]:
                if nbr not in deleted and v < nbr:
                    union(v, nbr)

        return len(set(find(v) for v in graph if v not in deleted))

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
        sd = [0] * n

        for v in range(n):
            if not G.adj_list[v]:
                continue
            graph, nVerts = kCoreBaseStructuralDiversity.subgraph(G, v)
            graph, nVerts, deleted = kCoreBaseStructuralDiversity.kCore(
                nVerts, graph, k)
            sd[v] = kCoreBaseStructuralDiversity.nKCores(graph, nVerts, deleted)

        return sd

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
