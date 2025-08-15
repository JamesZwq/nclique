#!/usr/bin/env python3
# Auto-generated for 5473837

STUDENT_ID = "5473837"
STUDENT_NAME = "Yuming Jiang"

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
        sd = [0] * n  # Initialize the output list for τ_k(v) values

        for v in range(n):
            neighbors = G.adj_list[v]  # Get neighbors of vertex v
            if not neighbors:
                continue

            # Step 1: Build the neighbor-induced subgraph
            induced_adj = kCoreBaseStructuralDiversity._induce_subgraph(G, neighbors)

            # Step 2: Extract the k-core from the induced subgraph
            k_core_adj = kCoreBaseStructuralDiversity._extract_k_core(induced_adj, k)

            # Step 3: Count the number of connected components in the k-core
            num_components = kCoreBaseStructuralDiversity._count_connected_components(k_core_adj)

            sd[v] = num_components  # τ_k(v) = number of connected k-core components

        return sd

    @staticmethod
    def _induce_subgraph(G, nodes):
        node_set = set(nodes)
        adj = defaultdict(set)
        for u in nodes:
            for v in G.adj_list[u]:  # Traverse edges only within the node set
                if v in node_set:
                    adj[u].add(v)
        return adj  # Return the induced subgraph as an adjacency list

    @staticmethod
    def _extract_k_core(adj, k):
        degree = {u: len(adj[u]) for u in adj}
        queue = deque([u for u in adj if degree[u] < k])

        # Remove nodes with degree less than k iteratively
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                adj[v].discard(u)
                if degree[v] >= k:
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)
            adj[u].clear()
            degree[u] = 0

        # Return only nodes that are part of the final k-core
        return {u: neighbors for u, neighbors in adj.items() if neighbors}

    @staticmethod
    def _count_connected_components(adj):
        visited = set()
        count = 0

        # Standard BFS to find connected components
        for node in adj:
            if node not in visited:
                queue = deque([node])
                visited.add(node)

                while queue:
                    u = queue.popleft()
                    for v in adj[u]:
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)
                count += 1

        return count



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
