#!/usr/bin/env python3
# Auto-generated for 5324147

STUDENT_ID = "5324147"
STUDENT_NAME = "(Evelyn) Meitong Zhou"

# ======= 学生代码 =======
from collections import deque

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
        neighbors = G.adj_list if hasattr(G, 'adj_list') else G.adjustList
        tau = [0 for _ in range(n)]  # Initialize result list

        def count_components(core_nodes, sub_adj):
            """Count the number of connected components in the induced subgraph."""
            vis = set()
            count = 0
            for u in core_nodes:
                if u not in vis:
                    count += 1
                    stack = [u]
                    vis.add(u)
                    while stack:
                        curr = stack.pop()
                        for w in sub_adj[curr]:
                            if w in core_nodes and w not in vis:
                                vis.add(w)
                                stack.append(w)
            return count

        for v in range(n):
            # Get neighbor list for vertex v (excluding self-loops)
            neighbors_v = [w for w in neighbors[v] if w != v]  # Exclude self-loops from neighbor list
            if len(neighbors_v) < k:
                tau[v] = 0
                continue
            # If not enough neighbors, k-core is empty for this node

            # Build neighbor-induced subgraph with unique undirected edges only, remove self-loops
            neighbors_set = set(neighbors_v)
            sub_adj = {node: [] for node in neighbors_set}
            for u in neighbors_set:
                for w in neighbors[u]:
                    if w in neighbors_set and u < w:   # only add once, prevent self-loop and duplicates
                        sub_adj[u].append(w)
                        sub_adj[w].append(u)
            # Construct the induced subgraph for v's neighbors (undirected, no self-loops, no duplicate edges)

            # K-core peeling
            candidates = set(neighbors_set)
            deg = {node: len(sub_adj[node]) for node in candidates}
            q = deque([node for node in candidates if deg[node] < k])
            removed = set()
            while q:
                u = q.popleft()
                if u in removed:
                    continue
                removed.add(u)
                for w in sub_adj[u]:
                    if w in candidates and w not in removed:
                        deg[w] -= 1
                        if deg[w] < k:
                            q.append(w)

            # Get nodes remain in k-core
            core_nodes = [node for node in candidates if node not in removed]
            if not core_nodes:
                tau[v] = 0
                continue

            # Count connected components in k-core
            tau[v] = count_components(core_nodes, sub_adj)
        return tau

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
