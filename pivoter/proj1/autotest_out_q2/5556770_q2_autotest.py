#!/usr/bin/env python3
# Auto-generated for 5556770

STUDENT_ID = "5556770"
STUDENT_NAME = "Xinyi Shen"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################
class kCoreBaseStructuralDiversity(object):
  def __init__(self):
    pass

  @staticmethod
  def dfs(subgraph, visited, start, node_set):
    stack = []  #used to simulate recursive DFS traversal
    stack.append(start)
    visited[start] = True
    while stack:
      current = stack.pop()
      for n in subgraph[current]:
        if n in node_set and not visited[n]:
          visited[n] = True
          stack.append(n)

  @staticmethod
  def build_subneighbour(G, v):
    "build the subgraph"
    neighbors = G.adj_list[v]
    id_map = dict()  #create a map from original vertex to new
    subgraph = []  #initialize the adjacency list for the neighbor-induced subgraph
    for new, ori in enumerate(neighbors):
        id_map[ori] = new
    for _ in range(len(neighbors)):
        subgraph.append([])
    for u in neighbors:
      for w in G.adj_list[u]:
        if w in id_map:
          subgraph[id_map[u]].append(id_map[w])
    return subgraph

  @staticmethod
  def peel_kcore(subgraph, k):
    "k-core decomposition"
    n = len(subgraph)
    degree = []  #initialize degree list for each node
    removed = [False] * n
    changed = True
    for l in subgraph:
        node_degree = len(l)
        degree.append(node_degree)
    while changed:
      changed = False
      for i in range(n):
        if not removed[i] and degree[i] < k:
          removed[i] = True
          changed = True
          for neighbor in subgraph[i]:
            degree[neighbor] -= 1
    return [e for e in range(n) if not removed[e]]

  @staticmethod
  def count_components(subgraph, nodes):
    "count connected components among remaining k-core nodes"
    count = 0
    n_set = set(nodes)
    length = len(subgraph)
    visited = [False] * length  #check if visited or not
    for n in n_set:
      if not visited[n]:
        kCoreBaseStructuralDiversity.dfs(subgraph, visited, n, n_set)
        count += 1
    return count

  @staticmethod
  def process(G, k):
    "compute the k-core-based structural diversity"
    n = G.vertex_num  #total number of vertices
    result = [0] * n
    for v in range(n):
      neighbors = G.adj_list[v]
      if not neighbors:
        continue
      subgraph = kCoreBaseStructuralDiversity.build_subneighbour(G, v)
      core_node = kCoreBaseStructuralDiversity.peel_kcore(subgraph, k)
      if core_node:
        result[v] = kCoreBaseStructuralDiversity.count_components(subgraph, core_node)
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
