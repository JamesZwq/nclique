#!/usr/bin/env python3
# Auto-generated for 5390453

STUDENT_ID = "5390453"
STUDENT_NAME = "Xiao Han"

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

      for v in range(n):  # Traverse each point v
        neighbors = G.adj_list[v]  # The neighbor list N(v) of the current point v does not contain v
        if not neighbors:  # If this point has no neighbors (is an isolated point)

          sd[v] = 0  # Then = 0, the structural diversity is 0
          continue

        nbr_graph = kCoreBaseStructuralDiversity.build_neighbor_subgraph(G, neighbors) #Construct an adjacency subgraph: only contain the points in neighbors and the edges between them

        sd[v] = kCoreBaseStructuralDiversity.compute_kcore_components(nbr_graph, k)


      return sd

    @staticmethod
    def build_neighbor_subgraph(G, neighbors):
      nbr_nodes = set(neighbors)                   # Use sets to remove duplicates and speed up searches (to facilitate subsequent determination of whether w is a neighbor)）
      nbr_graph = {u: [] for u in nbr_nodes}           # Initialize the adjacency subgraph: each neighbor point u corresponds to an empty adjacency list
      for u in nbr_nodes:                       # Traverse each neighbor node u (point in the subgraph)
        for w in G.adj_list[u]:                    # Traverse all neighbors w of u in the original graph
          if w in nbr_nodes:                     # If w is also a point in the subgraph (also in neighbors)
            nbr_graph[u].append(w)                   # Keep the edge u-w in the subgraph
      return nbr_graph                            # Returns the constructed adjacency subgraph nbr_graph (which is an adjacency list dictionary)


    @staticmethod
    def compute_kcore_components(nbr_graph, k):                # Input: adjacency subgraph, k value
      degrees = {}                              # Create an empty dictionary to hold the degree of each point
      for u in nbr_graph:                           # Traverse each point in the subgraph
        degrees[u] = len(nbr_graph[u])                     # Store the degree of point u (the number of its neighbors)

      queue = deque()                              # Create an empty skinning queue
      for u in degrees:                             # Iterate over the degree of each point
        if degrees[u] < k:                           # If its degree is less than k, it needs to be peeled off
          queue.append(u)                             # Put into the peeling queue


      while queue:                            # Start k-core peeling
        u = queue.popleft()                      # Take out a point that is not enough
        for v in nbr_graph[u]:                    # Traverse its neighbors
          if v in degrees:                     # The neighbor is still in the picture (not stripped off) This is a protection mechanism
            degrees[v] -= 1                    # The neighbor loses one connection, and the degree decreases by one
            if degrees[v] == k - 1:               # If the neighbor degree also becomes < k
              queue.append(v)                   # Join the skinning queue
        del degrees[u]                         # Remove u from the graph


      # At this point, the remaining points in degrees are the points in k-core

      visited = set()  # Used to record visited points
      count = 0  # Used to record the number of connected blocks

      for u in degrees:                            # Traverse each point that still exists in the k-core subgraph
        if u not in visited:                        # If it has not been visited yet, it means a new connected block is found.
          count += 1                            # Connected Block +1
          stack = [u]                            # Expand this connected block starting from u
          while stack:                            # As long as the stack is not empty, it will continue to expand
            curr = stack.pop()                      # The point currently being processed
            if curr in visited:                     # If you have already visited, skip
                continue
            visited.add(curr)                        # Otherwise mark as visited
            for nei in nbr_graph[curr]:                 # Traverse the neighbors of the current point
              if nei in degrees and nei not in visited:      # Neighbors are in the subgraph and have not been visited
                stack.append(nei)                    # Push into the stack and wait for later access

      return count  # This is τ_k(v) of the current adjacency subgraph: the number of connected blocks in the k-core





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
