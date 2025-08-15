#!/usr/bin/env python3
# Auto-generated for 5468765

STUDENT_ID = "5468765"
STUDENT_NAME = "Yuchen Bian"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################


class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _construct_induced_subgraph(G, v_neighbors):
        # construct G[N[v]], neighbor induced subgraph
        induced_subgraph = {}
        for neighbor in v_neighbors:
            induced_subgraph[neighbor] = set()

        # convert to set, faster than traverse list
        neighbor_set = set(v_neighbors)
        # for each neighbor
        for u in v_neighbors:
            # traverse this neighbor's neighbors
            for w in G.adj_list[u]:
                # if it in neighbor set and this is not a cycle
                if w in neighbor_set and w != u:
                    induced_subgraph[u].add(w)

        return induced_subgraph

    @staticmethod
    def _init_for_decomposition(induced_subgraph):
        # degree sorted array
        D = sorted(induced_subgraph.keys(),
                   key=lambda x: len(induced_subgraph[x]))

        # find max degree
        max_degree = 0
        if len(D) != 0:
            for node in D:
                node_degree = len(induced_subgraph[node])
                if node_degree > max_degree:
                    max_degree = node_degree
        # count number of nodes that each degree has
        bin_count = [0] * (max_degree + 2)
        for i in D:
            bin_count[len(induced_subgraph[i])] += 1
        # bin array, compute start position of each degree
        b = [0] * (max_degree + 2)
        for i in range(1, max_degree + 1):
            b[i] = b[i - 1] + bin_count[i - 1]

        # position array, store position of each node in D[]
        p = {}
        for i in range(len(D)):
            p[D[i]] = i

        # degree array, store current degree of each node
        d = {}
        for key, value in induced_subgraph.items():
            d[key] = len(value)

        return D, b, p, d

    @staticmethod
    def _swap(D, b, p, v, degree_v):
        # move v from degree[d[v]+1] to degree[d[v]]
        # v's position in D
        position_v = p[v]

        # start position of current degree
        position_w = None
        if degree_v < len(b):
            position_w = b[degree_v]
        else:
            position_w = len(D)

        w = D[position_w]

        if v != w:
            # swap v and w in D
            temp = D[position_v]
            D[position_v] = D[position_w]
            D[position_w] = temp

            temp = p[v]
            p[v] = p[w]
            p[w] = temp

        # update b
        b[degree_v] += 1

    @staticmethod
    def _core_decomposition(induced_subgraph, k):
        # init
        D, b, p, d = kCoreBaseStructuralDiversity._init_for_decomposition(
            induced_subgraph)

        if len(D) == 0:
            return []

        for i in range(len(D)):
            # already >= k, skip this node
            u = D[i]
            if d[u] >= k:
                continue

            # traverse u's neighbors
            for v in induced_subgraph[u]:
                # only update neighbors that degree is bigger than u
                if d[v] > d[u]:
                    degree_v = d[v]
                    # neighbor's degree - 1
                    d[v] -= 1
                    # swap
                    kCoreBaseStructuralDiversity._swap(D, b, p, v, degree_v)

        # store remain nodes, set is faster
        remain_nodes = set()
        for u in D:
            if d[u] >= k:
                remain_nodes.add(u)

        return remain_nodes

    @staticmethod
    def _bfs(remain_nodes, induced_subgraph):
        if len(remain_nodes) == 0:
            return 0

        # find the number of 't'
        result = 0
        visited = set()

        # bfs to each node
        for u in remain_nodes:
            # if this node is not visited, this is a start of new part
            if u not in visited:
                # init bfs
                queue = deque()
                visited.add(u)
                queue.append(u)
                # we find a new start, so this is a new part, update result
                result += 1

                # bfs
                while queue:
                    v = queue.popleft()
                    for w in induced_subgraph[v]:
                        if w in remain_nodes and w not in visited:
                            visited.add(w)
                            queue.append(w)

        return result

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

        # traverse all vertices
        for v in range(n):
            # find v's neighbors, N[v]
            v_neighbors = G.adj_list[v]

            # construct induced subgraph
            induced_subgraph = kCoreBaseStructuralDiversity._construct_induced_subgraph(
                G, v_neighbors)

            # k-core decomposition
            remain_nodes = kCoreBaseStructuralDiversity._core_decomposition(
                induced_subgraph, k)

            # find the number of connected subgraph
            result = kCoreBaseStructuralDiversity._bfs(
                remain_nodes, induced_subgraph)

            sd[v] = result

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
