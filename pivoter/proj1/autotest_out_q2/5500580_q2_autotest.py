#!/usr/bin/env python3
# Auto-generated for 5500580

STUDENT_ID = "5500580"
STUDENT_NAME = "Lu Zhao"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num  # Total number of vertices in the graph
        result = [0] * n  # Store the τ_k value of each vertex, initialized to 0

        # Traverse each node v in the graph and process them one by one
        for v in range(n):
          neighbors = G.adj_list[v]  # Get the neighbor list of node v
          if not neighbors:
              # If v is an isolated node (has no neighbors), then τ_k(v) = 0
              result[v] = 0
              continue

          # Step 1: Construct the neighbor subgraph nbr_v
          subgraph_nodes = set(neighbors)  # Adjacent point set
          subgraph_adj = defaultdict(set)  # Adjacency list of subgraph

          # Only keep edges between neighbors
          for u in subgraph_nodes:
              for w in G.adj_list[u]:
                  if w in subgraph_nodes:
                      subgraph_adj[u].add(w)

          # Step 2: Perform k-core deletion on the subgraph
          # Initialize the degree of each node
          degrees = {node: len(subgraph_adj[node]) for node in subgraph_nodes}
          # Add all nodes with a degree less than k to the queue and prepare for deletion
          queue = deque([node for node in subgraph_nodes if degrees[node] < k])

          # Start the k-core removal process
          while queue:
              u = queue.popleft()
              for nei in subgraph_adj[u]:
                  if degrees[nei] > 0:
                      degrees[nei] -= 1  # Update neighbor node degrees
                      if degrees[nei] == k - 1:
                          queue.append(nei)  # If the degree of a neighbor becomes k-1, it means that it no longer meets the k-core requirement
              degrees[u] = 0  # The current node has been removed and the degree is set to 0

          # The remaining nodes are valid nodes with degree ≥ k. Extract the remaining nodes
          valid_nodes = set(node for node in subgraph_nodes if degrees[node] >= k)

          if not valid_nodes:
              # All neighbor nodes are deleted and no k-core exists
              result[v] = 0
              continue

          # Step 3: Count the number of connected components consisting of valid nodes, i.e. the number of k-cores
          visited = set()
          count = 0

          # Use DFS to traverse unvisited connected components
          for node in valid_nodes:
              if node not in visited:
                  count += 1  # New connected component
                  # Use DFS to traverse the connected component
                  stack = [node]
                  while stack:
                      u = stack.pop()
                      if u in visited:
                          continue
                      visited.add(u)
                      # Traverse all neighbors connected to u that are still in the k-core
                      for nei in subgraph_adj[u]:
                          if nei in valid_nodes and nei not in visited:
                              stack.append(nei)

          # Save the number of connected components as the structural diversity value of the node
          result[v] = count

        return result


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
