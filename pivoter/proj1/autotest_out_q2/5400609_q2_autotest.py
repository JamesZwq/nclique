#!/usr/bin/env python3
# Auto-generated for 5400609

STUDENT_ID = "5400609"
STUDENT_NAME = "Haoyue Bu"

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

        for v in range(n):          # Traverse each vertex v.
            neighbours = G.adj_list[v]  # Get a list of v's neighbors (directly connected nodes).

            if len(neighbours) == 0:
                # If the vertex has no neighbors, it is an isolated point.
                continue

            set_neighbours = set(neighbours)  # Convert the neighbor list into a set to facilitate subsequent operations.

            if v in set_neighbours:
                # Prevent self-loops (v is connected to itself) and remove it proactively.
                set_neighbours.remove(v)

            if len(set_neighbours) == 0:
                # If the neighbor set is empty after removing the self-loop, it means that there is no neighbor subgraph to build, so skip.
                continue

            # Construct a neighbor-induced subgraph sub_G.
            sub_G = {}

            for u in set_neighbours:
                # Safety check to Skip unusual indexes.
                if u >= len(G.adj_list):
                    continue

                sub_G[u] = []  #Initialize the adjacency list of u.

                for w in G.adj_list[u]:
                    # Only add adjacent vertices in the neighbor set, ensuring that only connections within the neighbor subgraph are processed.
                    if w in set_neighbours and w != u:
                        sub_G[u].append(w)

            # Calculate the number of k-core connected blocks in the neighbor induced subgraph of the current vertex.
            sd[v] = kCoreBaseStructuralDiversity._count_k_core(sub_G, k)

        return sd


    @staticmethod
    def _count_k_core(sub_G, k):

        if not sub_G:
            return 0  # If the subgraph is empty, return 0 directly.

        # Initialize the degree (number of neighbors) of each node.
        degrees = {}
        for u in sub_G:
            degrees[u] = len(sub_G[u])

        deleted = set()      # Stores all nodes that will be deleted.
        queue = deque()

        #add all nodes with degree less than k to the deletion queue.
        for u in sub_G:
            if degrees[u] < k:
                queue.append(u)
                deleted.add(u)

        # k-core stripping: keep removing points with degree < k until stability.
        while queue:
            u = queue.popleft()
            for v in sub_G[u]:
                if v in deleted:
                    continue  # Skip deleted nodes
                degrees[v] -= 1  # Neighbor degree reduction
                if degrees[v] < k:
                    # Once degree < k, add to the deletion list.
                    deleted.add(v)
                    queue.append(v)

        # After stripping, the remaining nodes are the ones that truly belong to k-core.
        remains = []
        for u in sub_G:
            if u not in deleted:
                remains.append(u)

        if len(remains) == 0:
            # If there are no remaining nodes, there is no k-core structure.
            return 0

        visited = set()        # Mark visited nodes to prevent repeated traversal
        k_core_count = 0       # Count the number of k-core connected blocks

        # Perform BFS on the remaining nodes and count the number of connected blocks.
        for node in remains:
            if node in visited:
                continue  # Visited Skip

            # Start BFS from the starting point of a new connected block
            queue = deque([node])
            visited.add(node)
            k_core_count += 1  # Find a new connected block

            while queue:
                u = queue.popleft()
                for v in sub_G[u]:
                    # Only visit nodes that are still in k-core and have not been visited
                    if v not in visited and v not in deleted:
                        visited.add(v)
                        queue.append(v)

        return k_core_count  # Returns the final number of connected blocks

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
