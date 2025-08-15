#!/usr/bin/env python3
# Auto-generated for 5512062

STUDENT_ID = "5512062"
STUDENT_NAME = "Ziqi Zhou"

# ======= 学生代码 =======
from collections import deque
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        # Initialize the class
        pass
    @staticmethod
    def process(G, k):
        """
        Compute the k-core-based sd for each vertex in the graph.
        input:
        G : UndirectedUnweightedGraph
        k : int, the k value for k-core computation.
        output:
        List[int]  # τ_k(v) for all v
        """
        # Initialize the array with 0 for each vertex
        n = G.vertex_num
        sd = [0 for _ in range(n)]
        # Iterate through each vertex in the graph
        for v in range(n):
            all_neighbors = G.adj_list[v]
            # If the vertex has no neighbors, the sd is 0
            if not all_neighbors:
                sd[v] = 0
                continue
            # Construct the neighbor-induced subgraph (sub_G) containing only neighbors
            # and edges between them
            sub_G = {}
            for u in all_neighbors:
                sub_G[u] = []
            for u in all_neighbors:
                # Add edges to sub_G if both endpoints are in the neighbor set
                for a in G.adj_list[u]:
                    if a in all_neighbors:
                        sub_G[u].append(a)
            # Compute the number of k-cores
            sd[v] = kCoreBaseStructuralDiversity.count(sub_G, k)
        print(sd)
        return sd
    @staticmethod
    def count(sub_G, k):
        # If sub_G is empty, return 0
        if not sub_G:
            return 0
        # Initialize degree dictionary for each vertex in the subgraph
        degrees = {u: len(all_neighbors) for u, all_neighbors in sub_G.items()}
        # Initialize set to track vertices removing during k-core decomposition
        removed = set()
        # Initialize queue for vertices with degree less than k
        vertices = deque()
        for u in sub_G:
            if k > degrees[u]:
                vertices.append(u)
                removed.add(u)
        # Do k-core decomposition by removing vertices with degree < k
        while vertices:
            u = vertices.popleft()
            for v in sub_G[u]:
                if v not in removed:
                    degrees[v] = degrees[v] - 1
                    # If a neighbor's degree is smaller than k, mark it for removal
                    if k > degrees[v]:
                        removed.add(v)
                        vertices.append(v)

        # Identify remaining vertices (those in the k-core)
        remain_list = [u for u in sub_G if u not in removed]


        if not remain_list:
            return 0

        # Count k-cores in the subgraph
        visited = set()
        count1 = 0
        for node in remain_list:
            # using BFS
            if node not in visited:
                count1 = count1 +1
                queue = deque()
                queue.append(node)
                visited.add(node)
                # Explore the connected component
                while queue:
                    u = queue.popleft()
                    for v in sub_G[u]:
                        if v not in visited and v not in removed:
                            visited.add(v)
                            queue.append(v)
        return count1

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
