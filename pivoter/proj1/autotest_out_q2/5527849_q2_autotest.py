#!/usr/bin/env python3
# Auto-generated for 5527849

STUDENT_ID = "5527849"
STUDENT_NAME = "Shiyi Chen"

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
        n = G.vertex_num
        tau_k = [0] * n

        # do the global k_core decomposition
        core_numbers = kCoreBaseStructuralDiversity._k_core_decomposition(G)

        # only consider when core number >= k for later filters
        filter_res = [[] for _ in range(n)]
        for v in range(n):
            if core_numbers[v] >= k:
                filter_res[v] = [u for u in G.adj_list[v] if core_numbers[u] >= k]

        # use cahce to save the redundent calculations
        cache = {}
        # for every node, we compute the number of connected k_core number
        for v in range(n):
            if core_numbers[v] < k: # skip if core numbers less then k
                continue
            nbrs = filter_res[v]
            if not nbrs: # skip if empty
                continue

            key = frozenset(nbrs) # construct the key
            if key in cache: # if cache hit, use the result from the cache
                tau_k[v] = cache[key]
                continue

            # if not hit and not corner case, we construct the neighbour-induced subgraph
            nbrs_set = set(nbrs)
            induced_adj = defaultdict(list)
            for u in nbrs:
                induced_adj[u] = [w for w in filter_res[u] if w in nbrs_set]
            # compute the number of connected k_core in the subgraph
            connected_count = kCoreBaseStructuralDiversity._compute_k_core_num(nbrs, induced_adj, k)
            cache[key] = connected_count
            tau_k[v] = connected_count

        return tau_k

    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def _k_core_decomposition(G):
        n = G.vertex_num
        degrees = [len(G.adj_list[v]) for v in range(n)] # deg of each node
        core = degrees[:]
        max_deg = max(degrees) if n > 0 else 0

        # initialization
        class_deg = [[] for _ in range(max_deg + 1)]
        flag = [False] * n

        # for every node, add its degree to corrsponding class
        for v in range(n):
            class_deg[degrees[v]].append(v)
            flag[v] = True

        # k-core peeling when d less than max degree
        d = 0
        while d <= max_deg:
            while class_deg[d]:
                v = class_deg[d].pop()
                if not flag[v]:
                    continue
                flag[v] = False
                core[v] = d
                # update the degree of its neighbours
                for u in G.adj_list[v]:
                    if flag[u] and core[u] > d:
                        core[u] -= 1
                        class_deg[core[u]].append(u)
            d += 1
        return core

    @staticmethod
    def _compute_k_core_num(graph_nodes, graph_adj, k):
        if not graph_nodes:
            return 0
        # initialzation
        degrees = {node: len(graph_adj.get(node, [])) for node in graph_nodes}
        active_nodes = set(graph_nodes)
        # add node with degree less than k to the queue
        peel_nodes = deque(node for node, deg in degrees.items() if deg < k)

        # run BFS to do peeling process
        while peel_nodes:
            u = peel_nodes.popleft()
            if u not in active_nodes:
                continue
            active_nodes.remove(u)
            for v in graph_adj.get(u, []):
                if v in active_nodes and degrees[v] >= k:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        peel_nodes.append(v)

        # traverse every node in the k_core
        # if the node is not visted, a new connected k_core
        visited = set()
        num_k_cores = 0
        for node in active_nodes:
            if node not in visited:
                num_k_cores += 1
                track_nodes = [node]
                visited.add(node)
                # add neighbours from this connected k_core to be checked later
                while track_nodes:
                    cur = track_nodes.pop()
                    for nbrs in graph_adj.get(cur, []):
                        if nbrs in active_nodes and nbrs not in visited:
                            visited.add(nbrs)
                            track_nodes.append(nbrs)
        return num_k_cores

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
