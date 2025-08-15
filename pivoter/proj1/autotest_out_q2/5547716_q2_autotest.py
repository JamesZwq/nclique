#!/usr/bin/env python3
# Auto-generated for 5547716

STUDENT_ID = "5547716"
STUDENT_NAME = "Nguyen Phan"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # Total: O(n + m)
    @staticmethod
    def getKcore(G, active_vertices, k):
      # Init the degree
      # O(n)
      degree = [0] * G.vertex_num
      for u in range(G.vertex_num):
        for v in G.adj_list[u]:
          degree[u] += 1

      # Add the vertices with degree smaller than k
      # O(n)
      queue = deque()
      for u, d_u in enumerate(degree):
        if degree[u] < k:
          active_vertices[u] = False
          queue.append(u)

      # Iteratively remove nodes with degree < k
      # O(n + m)
      while len(queue) > 0:
        u = queue.popleft()
        for v in G.adj_list[u]:
          if not active_vertices[v]:
            continue
          degree[v] -= 1
          if degree[v] < k and active_vertices[v]:
            queue.append(v)
            active_vertices[v] = False

      return active_vertices

    @staticmethod
    def dfs(G, active_vertices, visited, s):
      if visited[s]:
        return

      stack = [s]
      while stack:
        u = stack.pop()
        visited[u] = True
        for v in G.adj_list[u]:
          if not active_vertices[v] or visited[v]:
            continue
          stack.append(v)


    # Total: O(n + m)
    @staticmethod
    def countKCore(G, active_vertices, k):
      # Get the vertices that are k-core
      # O(n + m)
      kCoreBaseStructuralDiversity.getKcore(G, active_vertices, k)

      # Count the number of groups in the remaining vertices
      # O(n + m)
      k_core_num = 0
      visited = [False] * G.vertex_num
      for u in range(G.vertex_num):
        if visited[u] or not active_vertices[u]:
          continue
        kCoreBaseStructuralDiversity.dfs(G, active_vertices, visited, u)
        k_core_num += 1

      return k_core_num

    # Total: O(n(n + m))
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
        # O(n)
        n = G.vertex_num
        sd = [0] * n

        # O(n) * O(n + m) = O(n(n + m))
        for u in range(n):
          # assigning new vertex index to u's neighbors
          # O(n)
          vertices_idx = [-1] * n;
          vertex_idx = 0
          for v in G.adj_list[u]:
            vertices_idx[v] = vertex_idx
            vertex_idx += 1

          # create the neighbor induced subgraph
          # O(n + m)
          edges_list = []
          active_vertices = []
          num_vertices = 0
          num_edges = 0
          for i in range(G.vertex_num):
              if vertices_idx[i] == -1:
                  continue
              active_vertices.append(True)
              num_vertices += 1
              for j in G.adj_list[i]:
                  if vertices_idx[j] == -1:
                      continue
                  if i < j:
                      edges_list.append([vertices_idx[i], vertices_idx[j]])
                      num_edges += 1
          edges_list.append([num_vertices, num_edges])
          edges_list[0], edges_list[-1] = edges_list[-1], edges_list[0]
          subgraph = UndirectedUnweightedGraph(edges_list)

          # get the number of k core from neighbor induced subgraph
          # O(n + m)
          sd[u] = kCoreBaseStructuralDiversity.countKCore(subgraph, active_vertices, k)

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
