#!/usr/bin/env python3
# Auto-generated for 5511934

STUDENT_ID = "5511934"
STUDENT_NAME = "Jiahao Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass
################################################################################
    # time complexity: m: edges; n: nodes; for node a, we have V' = deg(a), E' = ∑_{b∈N(a)} deg(b) (subgraph edges, worst)
    # _get_k_core: create degrees->O(V');traversal each edge->O(E') ==>O(V' + E')
    # _find_connected_components: traversal each node->O(V');BFS->O(V' + E') ==> O(V' + E')
    # _find_k_cores: O(V' + E') for each subgraph
    # process: for node a: O(deg(a) + V' + E' + V' + E') ==> O(E')
    # total: O(∑_{a∈V} E') = O(∑_{v∈V} deg(v)^2)
    #
################################################################################
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

        for a in range(n):
            nbrs = set(G.adj_list[a])
            if not nbrs:
                sd[a] = 0
                continue

            # create subgraph with neighbors of a
            subgraph = {}
            for b in nbrs:
                subgraph[b] = []
                for c in G.adj_list[b]:
                    if c in nbrs and c != b:
                        subgraph[b].append(c)

            # compute structural diversity for node a
            kcore_comps = kCoreBaseStructuralDiversity._find_k_cores(subgraph, k)
            sd[a] = len(kcore_comps)

        return sd

    @staticmethod
    def _find_k_cores(g, k):
        if not g:
            return []

        # get nodes in k-core
        kcore_nodes = kCoreBaseStructuralDiversity._get_k_core(g, k)
        if not kcore_nodes:
            return []

        # get connected components within k-core nodes
        comps = kCoreBaseStructuralDiversity._find_connected_components(g, kcore_nodes)
        return comps

    @staticmethod
    def _get_k_core(g, k):
        # initialize degrees and queue
        degrees = {x: len(g[x]) for x in g}
        queue = deque()
        removed = set()

        # enqueue nodes with degree < k
        for a in g:
            if degrees[a] < k:
                queue.append(a)
                removed.add(a)

        # remove low-degree nodes iteratively
        while queue:
            v = queue.popleft()
            for a in g[v]:
                if a not in removed:
                    degrees[a] -= 1
                    if degrees[a] < k:
                        queue.append(a)
                        removed.add(a)

        # remaining nodes form the k-core
        return set(g.keys()) - removed

    @staticmethod
    def _find_connected_components(g, nodes):
        if not nodes:
            return []

        visited = set()
        comps = []

        for a in nodes:
            if a not in visited:
                comp = set([a])
                queue = deque([a])
                visited.add(a)

                while queue:
                    b = queue.popleft()
                    for c in g[b]:
                        if c in nodes and c not in visited:
                            queue.append(c)
                            comp.add(c)
                            visited.add(c)

                comps.append(comp)
        return comps



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
