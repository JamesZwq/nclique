#!/usr/bin/env python3
# Auto-generated for 3399789

STUDENT_ID = "3399789"
STUDENT_NAME = "Weijia Zang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque,defaultdict
################################################################################



class kCoreBaseStructuralDiversity(object):

    @staticmethod
    def create_subgraph(neighbours, full_adj):
        """
        return a dict: global node -> its list of neighbors within the subgraph.
        """
        # remove duplicates to get the unique list of neighbors
        distinct_nbrs = set(neighbours)

        # adjacency list of the subgraph
        sub_adj = defaultdict(list)

        # outer loop iterates over each neighbor node in the local graph
        for neighbour in distinct_nbrs:
            # rnsure every node has an empty list
            sub_adj[neighbour]

            # retrieve the neighbor's neighbors from the global graph
            for nbr_of_nbrs in set(full_adj[neighbour]):
                # if the neighbor's neighbor is also within the local set,
                # it means the two are connected within the local scope
                if nbr_of_nbrs in distinct_nbrs:
                    # add the neighbor in the local adjacency list
                    sub_adj[neighbour].append(nbr_of_nbrs)

        # the local adjacency list now contains all connections in the induced subgraph
        return sub_adj


    @staticmethod
    def calculate_k_core_cc(sub_adj, k):
        """
        perform k-core peeling on sub_adj and return the number of connected components
        in the remaining graph.
        """
        # initialize degree of each vertex in the subgraph (remove duplicates using set)
        degree_dict = {}
        for vertex in sub_adj:
            degree_dict[vertex] = len(set(sub_adj[vertex]))

        # boolean dictionary to mark whether a node has been "peeled"
        deleted_dict = {}
        for vertex in sub_adj:
            deleted_dict[vertex] = False

        # initialize the peeling queue
        queue = deque()
        for vertex, degree in degree_dict.items():
            # enqueue all nodes with degree less than k for the first peeling round
            if degree < k:
                queue.append(vertex)

        # peeling process: dequeue nodes with degree < k and update neighbors' degrees
        while queue:
            vertex = queue.popleft()

            # skip if already peeled
            if deleted_dict[vertex]:
                continue

            # mark as peeled
            deleted_dict[vertex] = True

            # traverse all neighbors of the peeled node
            for neighbour in sub_adj[vertex]:
                # if neighbor is still active, reduce its degree by 1
                if not deleted_dict[neighbour]:
                    degree_dict[neighbour] -= 1
                    # if its degree just dropped from k to k-1, enqueue it for next round
                    if degree_dict[neighbour] == k - 1:
                        queue.append(neighbour)

        return kCoreBaseStructuralDiversity.calculate_components(sub_adj, deleted_dict)


    @staticmethod
    def calculate_components(sub_adj, deleted_dict):
        """
        count the number of connected components among the nodes remaining in k-core.
        """
        # after peeling, all remaining nodes are part of the k-core
        # perform BFS to count connected components
        visited = {}
        for vertex in sub_adj:
            visited[vertex] = False

        count = 0  # Counter for connected components

        # iterate over all nodes in the local graph
        for vertex in sub_adj:
            # if the node is still active and not visited, a new component is found
            if not deleted_dict[vertex] and not visited[vertex]:
                count += 1

                # start BFS from the current node, mark all connected active nodes as visited
                queue = deque([vertex])
                visited[vertex] = True

                while queue:
                    x = queue.popleft()
                    for y in sub_adj[x]:
                        if not deleted_dict[y] and not visited[y]:
                            visited[y] = True
                            queue.append(y)

        # after BFS finishes, count stores the number of connected components
        return count

    @staticmethod
    def process(G, k):
        """
        parameters:
        • G is the input undirected graph object with attributes:
          vertex_num (number of vertices) and adj_list (adjacency list).
        • k is the peeling threshold.
        • Return an integer list of length G.vertex_num.
        """

        n = G.vertex_num

        if k <= 0:
            return [0] * G.vertex_num

        full_adj = G.adj_list  # global adjacency list, reused for performance

        # number of k-core connected components in the induced subgraph of each vertex
        cc_number = [0] * n

        # Iterate through each vertex
        for v in range(n):
            neighbours = full_adj[v]  # the neighbor list of vertex v in the original graph

            # if no neighbors, the result is directly 0
            if not neighbours:
                cc_number[v] = 0
                continue

            # 1 build the induced subgraph from v’s neighbors,
            # only keeping edges between them
            sub_adj = kCoreBaseStructuralDiversity.create_subgraph(neighbours, full_adj)

            # 2 perform k-core peeling and count the number of connected components
            cc_number[v] = kCoreBaseStructuralDiversity.calculate_k_core_cc(sub_adj, k)

        return cc_number




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
