#!/usr/bin/env python3
# Auto-generated for 5542599

STUDENT_ID = "5542599"
STUDENT_NAME = "Shuaitian Liu"

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
        n = G.vertex_num  # Number of vertices in the graph
        sd = []  # Used to save the structural diversity value of each point (τ_k(v))

        # Process each point in the graph
        for v in range(n):
            nbr_set = set(G.adj_list[v])  # Get the neighbor set of the current point v

            # If the current point has no neighbors, then its structural diversity is 0
            if len(nbr_set) == 0:
                sd.append(0)
                continue  # Skip subsequent operations and process the next point

            # Step 1: Construct the edges of the neighbor-induced subgraph (consider only the edges between neighbors)
            induced_edges = set()  # Store the edges in the induced subgraph
            for u in nbr_set:
                neighbors_of_u = G.adj_list[u]  # u's neighbors
                for w in neighbors_of_u:
                    if w in nbr_set:  # If w is also a neighbor of v
                        if u < w:  # To avoid duplicate edges, only keep the smaller ones in front
                            edge = (u, w)
                            induced_edges.add(edge)

            # Step 2: Convert the edges of the induced subgraph into an adjacency list
            induced_adj = {}  # Adjacency list
            for node in nbr_set:
                induced_adj[node] = []  # Initialize the adjacency list for each neighbor
            for (u, w) in induced_edges:
                induced_adj[u].append(w)
                induced_adj[w].append(u)  # Because it is an undirected graph, both directions need to be added

            # Step 3: Execute k-core deletion operation
            remaining = set(nbr_set)  # Points currently in the graph
            changed = True  # Mark whether a deletion operation has occurred
            while changed:
                changed = False
                to_remove = set()  # Points to be deleted in this round
                for u in remaining:
                    degree = 0
                    for neighbor in induced_adj[u]:
                        if neighbor in remaining:  # Only count neighbors that have not been deleted
                            degree += 1
                    if degree < k:  # If the degree is less than k, prepare to delete
                        to_remove.add(u)
                        changed = True
                for u in to_remove:
                    remaining.remove(u)  # Remove these points from the graph

            # Step 4: Count the number of connected blocks of the remaining points (using BFS)
            visited = set()
            count = 0  # Number of connected blocks
            for node in remaining:
                if node not in visited:
                    count += 1  # New connected blocks
                    queue = deque()
                    queue.append(node)
                    visited.add(node)
                    while len(queue) > 0:
                        curr = queue.popleft()
                        for neighbor in induced_adj[curr]:
                            if neighbor in remaining:
                                if neighbor not in visited:
                                    visited.add(neighbor)
                                    queue.append(neighbor)

            # Add the τ_k value of the current point v to the result list
            sd.append(count)

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
