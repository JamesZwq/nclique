#!/usr/bin/env python3
# Auto-generated for 5501378

STUDENT_ID = "5501378"
STUDENT_NAME = "Jingyuan Wang"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core based structural diversity for all vertices in G.

        Parameters:
        G : UndirectedUnweightedGraph - the input graph
        k : int - the core number

        Returns:
        List[int] - τ_k(v) for all vertices v in G
        """
        n = G.vertex_num
        tau = [0] * n

        for v in range(n):
            # Get neighbors of v
            neighbors = G.adj_list[v]
            if not neighbors:
                tau[v] = 0
                continue

            # Create neighbor-induced subgraph
            neighbor_set = set(neighbors)
            subgraph_nodes = neighbors
            subgraph_adj = {}

            # Build adjacency list for subgraph
            for u in subgraph_nodes:
                subgraph_adj[u] = [nbr for nbr in G.adj_list[u] if nbr in neighbor_set]

            # Compute k-cores in the subgraph
            tau[v] = kCoreBaseStructuralDiversity.count_k_cores(subgraph_adj, k)

        return tau

    @staticmethod
    def count_k_cores(adj, k):
        """
        Count the number of k-cores in a subgraph.

        Parameters:
        adj : dict - adjacency list of the subgraph
        k : int - core number

        Returns:
        int - number of k-cores
        """
        if not adj:
            return 0

        # Make a copy of the adjacency list to work with
        degrees = {u: len(neighbors) for u, neighbors in adj.items()}
        remaining_nodes = set(adj.keys())

        # Initialize queue with nodes having degree < k
        queue = deque([u for u in remaining_nodes if degrees[u] < k])

        # Core decomposition
        while queue:
            u = queue.popleft()
            if u not in remaining_nodes:
                continue

            remaining_nodes.remove(u)

            for v in adj[u]:
                if v in remaining_nodes:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)

        # Now find connected components in remaining_nodes (these are k-cores)
        visited = set()
        core_count = 0

        for node in remaining_nodes:
            if node not in visited:
                # BFS to find connected component
                queue = deque([node])
                visited.add(node)
                core_count += 1

                while queue:
                    u = queue.popleft()
                    for v in adj[u]:
                        if v in remaining_nodes and v not in visited:
                            visited.add(v)
                            queue.append(v)

        return core_count

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
