#!/usr/bin/env python3
# Auto-generated for 5470334

STUDENT_ID = "5470334"
STUDENT_NAME = "Alexandra Ukrainskaya"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def count_components(nodes, edges):
      if len(nodes) == 0:
        return 0

      unvisited = set(nodes)
      count = 0

      while unvisited:
        count += 1
        start = unvisited.pop()
        queue = deque([start])

        while queue:
            node = queue.popleft()
            for neigh in edges[node]:
                if neigh in unvisited:
                    unvisited.remove(neigh)
                    queue.append(neigh)

      return count




    @staticmethod
    def count_k_cores(nodes, edges, k):
      #if k == 0:
      #  return kCoreBaseStructuralDiversity.count_zero_k_cores(nodes, edges)
      degrees = {node:len(edges[node]) for node in nodes}
      queue = deque([node for node in nodes if degrees[node] < k])
      is_visited = {node: degrees[node] < k for node in nodes}
      while queue:
        node = queue.popleft()
        for neigh in edges[node]:
          degrees[neigh] -= 1
          edges[neigh].remove(node)
          if degrees[neigh] < k:
            if not is_visited[neigh]:
              is_visited[neigh] = True
              queue.append(neigh)
        nodes.remove(node)
        edges.pop(node)
      return kCoreBaseStructuralDiversity.count_components(nodes, edges)


    @staticmethod
    def neighbour_induced_subgraph(adj_list, edge_list):
      nodes = set(adj_list) #used to be list
      edges = {}
      for node in adj_list:
        edges[node] = []

      for edge in edge_list:
        if edge[0] in nodes and edge[1] in nodes:
          edges[edge[0]].append(edge[1])
          edges[edge[1]].append(edge[0])

      return nodes, edges


    @staticmethod
    def get_edge_list(G):
      edge_list = []
      for v in range(G.vertex_num):
        for neigh in G.adj_list[v]:
          if neigh > v:
            edge_list.append((v, neigh))  #neigh, v

      return edge_list


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
      n = G.vertex_num
      neighbor_count = [len(G.adj_list[node]) for node in range(G.vertex_num)]
      sd = [0] * n
      #get edge list for faster computation of neighbour-induced subgraphs
      edge_list = kCoreBaseStructuralDiversity.get_edge_list(G)
      #create a dictionary for neighbout-induced subgraph
      for v in range(G.vertex_num):
        #get the subgraph
        if neighbor_count[v] >= k:
          subgraph_nodes, subgraph_edges = kCoreBaseStructuralDiversity.neighbour_induced_subgraph(G.adj_list[v], edge_list)
          num_k_cores = kCoreBaseStructuralDiversity.count_k_cores(subgraph_nodes, subgraph_edges, k)
          sd[v] = num_k_cores
        else:
          sd[v] = 0

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
