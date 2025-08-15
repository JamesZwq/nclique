#!/usr/bin/env python3
# Auto-generated for 5452789

STUDENT_ID = "5452789"
STUDENT_NAME = "Dan Jiang"

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
        tau = [0] * n

        # Process each vertex
        for v in range(n):
          neighbors = G.adj_list[v]

          if len(neighbors) == 0:
            tau[v] = 0
            continue

          # Build neighbor-induced subgraph
          neighbor_set = set(neighbors)
          vertex_map = {neighbors[i]: i for i in range(len(neighbors))}

          # Compute initial degrees in neighbor-induced subgraph
          degrees = [0] * len(neighbors)
          for i, u in enumerate(neighbors):
              for w in G.adj_list[u]:
                  if w in neighbor_set:
                      degrees[i] += 1

          # k-core decomposition: iteratively remove vertices with degree < k
          active = [True] * len(neighbors)
          queue = deque()

          # Initialize queue with vertices having degree < k
          for i in range(len(neighbors)):
            if degrees[i] < k:
                queue.append(i)
                active[i] = False

          # Process removals
          while queue:
              u_idx = queue.popleft()
              u = neighbors[u_idx]

              # Update degrees of neighbors of removed vertex
              for w in G.adj_list[u]:
                  if w in vertex_map and active[vertex_map[w]]:
                      w_idx = vertex_map[w]
                      degrees[w_idx] -= 1
                      if degrees[w_idx] < k:
                          active[w_idx] = False
                          queue.append(w_idx)

          # Count connected components in remaining graph
          visited = [False] * len(neighbors)
          components = 0

          for i in range(len(neighbors)):
              if active[i] and not visited[i]:
                  components += 1
                  # BFS to mark all vertices in this component
                  bfs_queue = deque([i])
                  visited[i] = True

                  while bfs_queue:
                      curr_idx = bfs_queue.popleft()
                      curr = neighbors[curr_idx]

                      # Check all neighbors of current vertex
                      for w in G.adj_list[curr]:
                          if w in vertex_map:
                              w_idx = vertex_map[w]
                              if active[w_idx] and not visited[w_idx]:
                                  visited[w_idx] = True
                                  bfs_queue.append(w_idx)

          tau[v] = components

        return tau


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
