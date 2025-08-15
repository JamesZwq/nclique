#!/usr/bin/env python3
# Auto-generated for 5526032

STUDENT_ID = "5526032"
STUDENT_NAME = "Shuaixian Tian"

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

        n: g.vartex_num, the total number of nodes in the graph
        m: Total number of edges in the graph
        d: Maximum degree of a single node
        """
        # TODO
        n = G.vertex_num
        sd = [0] * n
        # Iterate over all vertices to calculate the number of k-cores in the neighbour-induced subgraph for each vertex separately
        for node in range(n):
            neighbors = G.adj_list[node]
            if not neighbors:
                sd[node] = 0  # If the node has no neighbors, then the neighbor induced subgraph is empty
                continue
            neighbors_nodes = set(neighbors)  # Neighbor induced subgraph

            # compute τ_k(v) for current node v
            sd[node] = kCoreBaseStructuralDiversity.k_Core_function(G, neighbors_nodes, k)
            # The total time complexity is O(n*d^2) in the extreme case, and close to O(m) in the extremely simple and sparse case.

        return sd

    @staticmethod
    def k_Core_function(G, nodes, k):
        """
        Return: τ_k(v), the number of connected components in the k-core

        """
        if not nodes:
            return 0  # τ_k(v) = 0 if this node has no neighbors

        # Build neighbor induced subgraph (keep only connections inside nodes)
        subgraph = {}  # Initialize the empty dictionary
        for u in nodes:
            subgraph[u] = []
        for u in nodes:
            for v in G.adj_list[u]:
                if v in nodes:
                    subgraph[u].append(v) # Adjacency List structure: Only keep edges between neighbors
                    # Time complexity: outer d nodes, inner d neighbors, total O(d^2)


        # extract k-core (remove nodes with degree < k)
        degrees = {}  # Initialize the empty dictionary
        for u in subgraph:
            degrees[u] = len(subgraph[u])# calculate the degree of the current set of nodes

        queue = deque()  # Create an empty queue with the set of nodes to be removed
        for u in degrees:
            if degrees[u] < k:
                queue.append(u) # set of nodes to remove

        delete = set(queue)  # The set of deleted nodes

        while queue: # If the queue ends empty, it can't be deleted
            u = queue.popleft() # Extract from the left
            for v in subgraph[u]:
                if v not in delete:
                    degrees[v] -= 1 # update neighbor degrees
                    if degrees[v] < k:
                        delete.add(v) # update the set of removed nodes
                        queue.append(v) # update the set of nodes to remove
                        #Time complexity: Worst case every point and edge needs to be traversed, total O(d^2)


        # Get the remaining nodes in k-core
        remain = []  # Initialize the list of remaining nodes, subgraphs
        for u in subgraph:
            if u not in delete:
                remain.append(u)
                # Time complexity: Scan a node set O(d)

        if not remain:
            return 0  # τ_k(v) = 0 if there are no nodes left

        # Count number of connected components (BFS traversal)
        visited = set()
        count = 0  # Number of connected components

        for u in remain:
            if u not in visited:
                count += 1  # new connected block
                node_queue = deque([u])
                visited.add(u)

                while node_queue:
                    current_node = node_queue.popleft()# Extract from the left
                    for neighbor in subgraph[current_node]:
                        # Iterate over only undeleted and unvisited neighbors
                        if neighbor not in delete and neighbor not in visited:
                            visited.add(neighbor)
                            node_queue.append(neighbor)
        # A subgraph is connected if the entire subgraph is accessed in one BFS. Otherwise, there exist multiple connected graphs
        # Time complexity: If the deleted image is close to the complete graph, the time complexity can be close to O(d^2)
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
