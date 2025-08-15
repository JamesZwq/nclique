#!/usr/bin/env python3
# Auto-generated for 5510340

STUDENT_ID = "5510340"
STUDENT_NAME = "Jingqi Du"

# ======= 学生代码 =======
from collections import deque, defaultdict

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    ###########################################################################
    #Compute k-core-based structural diversity for each node in the graph
    def process(G, k):
        n = G.vertex_num
        result = [0] * n

        for v in range(n):
            adj_v = G.adj_list[v]
            if len(adj_v) < k:
                continue  # node degree less than k

            # Collect neighbors as a set
            neighbor_set = set(G.adj_list[v])
            induced = defaultdict(set)

            for node in G.adj_list[v]:
                for adj in G.adj_list[node]:
                    if adj in neighbor_set and node < adj:
                        induced[node].add(adj)
                        induced[adj].add(node)

            node_count = len(induced)
            if node_count < k:
                continue  # Skip if induced subgraph too small

            diversity = kCoreBaseStructuralDiversity._count_k_core(induced, k)
            result[v] = diversity

        return result
    @staticmethod
    @staticmethod
    ###########################################################################
    #Compute the number of connected components in the k-core of a given undirected subgraph
    def _count_k_core(subgraph, k):
        if len(subgraph) == 0:
            return 0

        # Initialize degree map and deletion set
        degree_map = {}
        #Build a mapping of node
        for node in subgraph:
            degree_map[node] = len(subgraph[node])
        to_remove = deque()
        #Initialize a queue with nodes
        for node in degree_map:
            if degree_map[node] < k:
                to_remove.append(node)
        removed = set()

        # Prune nodes with degree < k
        while to_remove:
          # Take one node from the queue
            target = to_remove.popleft()
            if target in removed:
                continue
            # Mark this node as removed from the k-core
            removed.add(target)
            neighbor_list = subgraph.get(target, set())
            # Traverse through each neighbor of the removed node
            for n in neighbor_list:
                if n not in removed:
                    new_deg = degree_map[n] - 1
                    degree_map[n] = new_deg
                    if new_deg < k:
                        to_remove.append(n)


        # Remaining nodes form the k-core
        core_nodes = list(filter(lambda x: x not in removed, subgraph))
        if not core_nodes:
            return 0

        # Count connected components in k-core
        visited = set()
        components = 0
        # Traverse each remaining node in the k-core
        for start in core_nodes:
            if start in visited:
                continue
            components += 1
            stack = [start]
            visited.add(start)
            # Perform DFS to visit all nodes in this connected component
            while stack:
                current = stack.pop()
                for neighbor in subgraph[current]:
                    if neighbor not in visited and neighbor not in removed:
                        visited.add(neighbor)
                        stack.append(neighbor)

        return components
# Return the total number of connected components in the k-core

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
