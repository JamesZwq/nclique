#!/usr/bin/env python3
# Auto-generated for 5523688

STUDENT_ID = "5523688"
STUDENT_NAME = "Yixue Yang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        for v in range(n):
            # Get the neighbor set of vertex v
            neighbors = set(G.adj_list[v])

            if len(neighbors) == 0:
                # If there are no neighbors, the structural diversity is 0
                sd[v] = 0
                continue

            # Construct neighbor induced subgraph
            # Neighbor list
            neighbor_nodes = list(neighbors)
            # Mapping of neighbor nodes to their indexes
            node_index_map = {}
            for idx, neighbor in enumerate(neighbor_nodes):
                node_index_map[neighbor] = idx

            # Construct the adjacency list of the neighbor-induced subgraph
            induced_adj = []
            n = len(neighbor_nodes)
            for _ in range(n):
                induced_adj.append([])
            for i, node in enumerate(neighbor_nodes):
                for adjacent in G.adj_list[node]:
                    if adjacent in node_index_map:
                        induced_adj[i].append(node_index_map[adjacent])

            # Find the number of all k-cores in the neighbor-induced subgraph
            sd[v] = kCoreBaseStructuralDiversity.count_k_cores(induced_adj, k)

        return sd

    @staticmethod
    def count_k_cores(graph_adj, k):
        """
        Calculate the number of k-cores in the graph.
        Use iteration to delete vertices with degree less than k, and then calculate the number of remaining connected components.
        """
        n = len(graph_adj)
        if n == 0:
            return 0

        # Calculate the degree of each vertex
        degrees = []
        for neighbors in graph_adj:
            degrees.append(len(neighbors))

        # Iterate and delete vertices with degree less than k
        remaining = set(range(n))
        queue = deque()

        # Initialization (add all vertices with degree less than k to the queue)
        for v in range(n):
            if degrees[v] < k:
                queue.append(v)

       # Iterate and delete vertices with degree less than k
        while queue:
            v = queue.popleft()
            if v not in remaining:
                continue

            remaining.remove(v)

            # Update neighbor degrees
            for neighbor in graph_adj[v]:
                if neighbor in remaining:
                    degrees[neighbor] = degrees[neighbor] - 1
                    if degrees[neighbor] < k:
                        queue.append(neighbor)

        if not remaining:
            return 0

        # Calculate the number of connected components directly in the remaining vertices
        visited = set()
        num_components = 0

        for v in remaining:
            if v not in visited:
                # BFS traverses connected components
                search_queue = deque([v])
                visited.add(v)

                while search_queue:
                    current_node = search_queue.popleft()
                    for neighbor in graph_adj[current_node]:
                        if neighbor in remaining:
                            if neighbor not in visited:
                                visited.add(neighbor)
                                search_queue.append(neighbor)

                num_components = num_components + 1

        return num_components

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
