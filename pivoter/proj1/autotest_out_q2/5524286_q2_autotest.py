#!/usr/bin/env python3
# Auto-generated for 5524286

STUDENT_ID = "5524286"
STUDENT_NAME = "Fei Ye"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # Compute the k-core subgraph of an undirected graph by removing nodes with degree < k until all remaining nodes have degree ≥ k.
    # Time Complexity: O(|V| + |E|)
    @staticmethod
    def compute_k_core(adj, k):
      # Calculate the current degree of each node
      # Time Complexity: O(|V| + |E|)
      degree = {v: len(neighbors) for v, neighbors in adj.items()}
      # Record whether the node is deleted
      # Time Complexity: O(|V|)
      deleted = {v: False for v in adj.keys()}
      # Initialize the queue: put all nodes with degree < k into the queue
      # Time Complexity: O(|V|)
      q = deque([v for v in adj.keys() if degree[v] < k])
      # Peeling process: continuously delete nodes with degree less than k
      # Perform at most |V| dequeues
      while q:
        # O(1)
        v = q.popleft()
        if deleted[v]:
          continue
        # Delete node v
        deleted[v] = True
        # For each neighbor u of v, update the degree and let u with degree less than k enter the queue
        # Iterate over all neighbors of v, ∑vdeg(v) = O(|E|)
        for u in adj[v]:
          if not deleted[u]:
            degree[u] = degree[u] - 1
            adj[u].remove(v)
            if degree[u] < k:
              q.append(u)
        # Clearing v’s neighbor list is a total of ∑vdeg(v) = O(|E|) operations
        adj[v].clear()
      # Build and return a k-core subgraph: traverse all nodes and their neighbors
      # Time complexity: O(|V| + |E|)
      k_core = {}
      for v in adj.keys():
        if not deleted[v]:
          neighbors = []
          for u in adj[v]:
            if not deleted[u]:
              neighbors.append(u)
          k_core[v] = neighbors

      return k_core

    # Count the number of connected components of an undirected graph
    # Time Complexity: O(|V| + |E|)
    @staticmethod
    def count_cc(adj):
      # Get all node sets
      # Time complexity: O(|V|) Need to traverse all vertices
      nodes = set(adj.keys())
      visited = set()
      count = 0
      # Loop: O(|V|) times
      for node in nodes:
        if node not in visited:
          queue = deque([node]) # O(1)
          visited.add(node) # O(1)
          while queue:
            v = queue.popleft() # O(1)
            # Iterate over all neighbors of v
            # Traverse all edges once cumulatively: ∑ deg(v) = O(|E|)
            for u in adj[v]:
              if u not in visited:
                visited.add(u)
                queue.append(u)
          count += 1
      return count

    # Time Complexity: O(|V|^3)
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
        # Construct the original graph adjacency list adj
        # Traverse all nodes and traverse their neighbors
        # total number of operations = ∑u deg(u) = O(|E|)
        adj = {i: [] for i in range(G.vertex_num)}
        for u in range(G.vertex_num):
            for v in G.adj_list[u]:
              adj[u].append(v)
        # Traverse each vertex and construct a neighbor induced subgraph for each vertex
        # Outer loop O(|V|) times
        for v in range(n):
          N_v = adj[v]
          if not adj[v]:
            sd[v] = 0
            continue
          # Construct the neighbor induced subgraph of v
          # O(|V|)
          Nset = set(N_v)
          nbr_adj = {u: set() for u in N_v}
          # Double traversal: for each neighbor u of v, traverse all neighbors w of u
          # Time Complexity: O(|V|^2)
          for u in N_v:
            for w in adj[u]:
              if w in Nset and w != u:
                  nbr_adj[u].add(w)
                  nbr_adj[w].add(u)
          # O(|V| + |E|)
          k_core = kCoreBaseStructuralDiversity.compute_k_core(nbr_adj, k)
          # O(|V| + |E|)
          num_cc = kCoreBaseStructuralDiversity.count_cc(k_core)
          sd[v] = num_cc
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
