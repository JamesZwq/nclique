#!/usr/bin/env python3
# Auto-generated for 5605054

STUDENT_ID = "5605054"
STUDENT_NAME = "Zhihao Wu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core-based structural diversity for all vertices

        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        # TODO
        n = G.vertex_num
        result = []

        # process each vertex
        for v in range(n):
            # get neighbors of v
            neighbors = G.adj_list[v]

            # if no neighbors, tau is 0
            if len(neighbors) == 0:
                result.append(0)
                continue

            # build subgraph of neighbors
            # use set for faster lookup
            neighbor_set = set(neighbors)

            # map neighbors to indices
            neighbor_list = list(neighbors)
            neighbor_map = {}
            for i in range(len(neighbor_list)):
                neighbor_map[neighbor_list[i]] = i

            # create adjacency list for subgraph
            subgraph_n = len(neighbor_list)
            subgraph_adj = []
            for i in range(subgraph_n):
                subgraph_adj.append([])

            # add edges between neighbors
            for u in neighbor_list:
                u_idx = neighbor_map[u]
                # check u's neighbors
                for w in G.adj_list[u]:
                    # only add if w is also neighbor of v
                    if w in neighbor_set and w != u:
                        w_idx = neighbor_map[w]
                        subgraph_adj[u_idx].append(w_idx)

            # count k-cores
            tau = find_k_cores(subgraph_adj, k)
            result.append(tau)

        return result


def find_k_cores(adj_list, k):
    """
    Find number of k-cores in graph
    """
    n = len(adj_list)

    if n == 0:
        return 0

    # special case k=0
    if k <= 0:
        return count_components(adj_list)

    # calculate degrees
    degrees = []
    active = []
    for i in range(n):
        degrees.append(len(adj_list[i]))
        active.append(True)

    # remove vertices with degree < k
    queue = deque()

    # add initial vertices
    for i in range(n):
        if degrees[i] < k:
            queue.append(i)

    # process queue
    while len(queue) > 0:
        u = queue.popleft()

        if active[u]:
            active[u] = False

            # update neighbors
            for v in adj_list[u]:
                if active[v]:
                    degrees[v] = degrees[v] - 1
                    if degrees[v] < k:
                        queue.append(v)

    # count components in remaining graph
    visited = [False] * n
    num_components = 0

    for i in range(n):
        if active[i] and not visited[i]:
            # DFS
            stack = [i]
            while stack:
                u = stack.pop()
                if not visited[u]:
                    visited[u] = True
                    for v in adj_list[u]:
                        if active[v] and not visited[v]:
                            stack.append(v)

            num_components = num_components + 1

    return num_components


def count_components(adj_list):
    """
    Count connected components
    """
    n = len(adj_list)
    if n == 0:
        return 0

    visited = [False] * n
    count = 0

    for i in range(n):
        if not visited[i]:
            # DFS from vertex i
            stack = [i]
            while stack:
                v = stack.pop()
                if not visited[v]:
                    visited[v] = True
                    for u in adj_list[v]:
                        if not visited[u]:
                            stack.append(u)
            count = count + 1

    return count

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
