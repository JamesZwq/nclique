#!/usr/bin/env python3
# Auto-generated for 5618743

STUDENT_ID = "5618743"
STUDENT_NAME = "Ruyun Xue"

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

        for v in range(n):
          subG = kCoreBaseStructuralDiversity.neighbour_induced_subgraph(G, v)

          # Step 1 & 2: Perform k-core decomposition on the neighbor-induced subgraph and directly get valid vertices
          valid = kCoreBaseStructuralDiversity.k_core_decomposition(subG, k)

          if not valid:
            sd[v] = 0
          else:
            sd[v] = kCoreBaseStructuralDiversity.count_components(subG, valid)

          # Step 3: Compute number of connected components on the valid subgraph

          # print(f"Number of connected components for vertex {v}: {component_count}")
          # sd[v] = component_count

        """
        if k == -1:  # Debug entry
          subG = kCoreBaseStructuralDiversity.neighbour_induced_subgraph(G, 2)
          print("Edge list (neighbor subgraph of vertex 2):")
          print(subG.edge_list)

          print("Adjacency list:")
          for i, nbrs in enumerate(subG.adj_list):
              print(f"{i}: {nbrs}")

          print("Debug finished, not returning τ_k array")
          return None  # or return []
        """
        return sd

    @staticmethod
    def neighbour_induced_subgraph(G, v):
      neighbours = set(G.adj_list[v])
      # Does not include v itself, only includes neighbors of v

      mapping = {old_id: new_id for new_id, old_id in enumerate(neighbours)}
      edges = []
      for u in neighbours:
        for w in G.adj_list[u]:
          if w in neighbours and mapping[u] < mapping[w]:
            edges.append([mapping[u], mapping[w]])
      vertex_num = len(mapping)
      edge_num = len(edges)
      edge_list = [[vertex_num, edge_num]] + edges

      """
      # Print debug info
      print("Neighbor set:", neighbours)
      print("Mapping:", mapping)
      print("Edge list:", edges)
      print("edge_list (for constructing graph):", edge_list)
      """
      return UndirectedUnweightedGraph(edge_list)

    @staticmethod
    def k_core_decomposition(G, k):
      n = G.vertex_num
      deg = [len(G.adj_list[v]) for v in range(n)]
      removed = [False] * n

      # Repeatedly remove vertices with degree less than k
      changed = True
      while changed:
        changed = False
        for v in range(n):
          if not removed[v] and deg[v] < k:
            removed[v] = True
            changed = True
            # Update the degrees of neighbors
            for u in G.adj_list[v]:
              if not removed[u]:
                deg[u] -= 1

      # Return remaining vertices (vertices in the k-core)
      k_core_vertices = [v for v in range(n) if not removed[v]]
      return k_core_vertices

    @staticmethod
    def count_components(G, valid):
        visited = [False] * G.vertex_num
        valid_set = set(valid)
        count = 0

        def dfs(v):
            visited[v] = True
            for u in G.adj_list[v]:
                if u in valid_set and not visited[u]:
                    dfs(u)

        for v in valid:
            if not visited[v]:
                dfs(v)  # Traverse the entire connected component
                count += 1  # Only increment after completing one component

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
