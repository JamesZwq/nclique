#!/usr/bin/env python3
# Auto-generated for 5493057

STUDENT_ID = "5493057"
STUDENT_NAME = "Zhicheng Han"

# ======= 学生代码 =======
class kCoreBaseStructuralDiversity(object):

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
    # Total number of vertices (assumed labeled 0,1,...,n-1)
    n = G.vertex_num

    # Bind the neighbor-listing function to nbr_func for clarity
    # Here we assume the graph implementation stores adjacency in G.adj_list
    nbr_func = lambda v: G.adj_list[v]

    # Initialize the result list with zeros
    sd = [0] * n

    # Process each vertex v in turn
    for v in range(n):
      # 1.Extract the list of neighbors of v
      nbrs = nbr_func(v)

      # 2.Build the adjacency structure for the induced subgraph on nbrs
      #  Represented as a dict: adj[u] = set of neighbors of u within nbrs
      adj = {u: set() for u in nbrs}
      for u in nbrs:
        # For each neighbor u of v, look at u's own neighbors
        for w in nbr_func(u):
          # Only keep edges between neighbors of v
          if w in adj:
            adj[u].add(w)

      # 3.Perform k-core pruning on this induced subgraph:
      #  Start with all nodes in core_nodes
      #  Compute initial degrees inside the induced subgraph
      core_nodes = set(adj.keys())
      deg = {u: len(adj[u]) for u in core_nodes}

      # Initialize a queue with any nodes whose degree < k
      queue = [u for u in core_nodes if deg[u] < k]

      # Iteratively remove nodes of degree < k, and decrement the degree of their neighbors (push neighbors that drop below k onto the queue)
      while queue:
        u0 = queue.pop()
        if u0 in core_nodes:
          core_nodes.remove(u0)
          for w in adj[u0]:
            if w in core_nodes:
              deg[w] -= 1
              if deg[w] == k - 1:
                queue.append(w)

      # 4.Count the number of connected components in the remaining k-core
      visited = set()
      comps = 0
      for u in core_nodes:
        if u not in visited:
          comps += 1
          # DFS search to mark all nodes in this component
          stack = [u]
          visited.add(u)
          while stack:
            x = stack.pop()
            for w in adj[x]:
              if w in core_nodes and w not in visited:
                visited.add(w)
                stack.append(w)

      # 5. Record the component count as τ_k(v)
      sd[v] = comps

    return sd

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
