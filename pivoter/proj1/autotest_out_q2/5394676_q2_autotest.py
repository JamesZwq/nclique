#!/usr/bin/env python3
# Auto-generated for 5394676

STUDENT_ID = "5394676"
STUDENT_NAME = "(Rachel) Xingyi Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################


class kCoreBaseStructuralDiversity(object):

    def __init__(self):
        pass

    @staticmethod
    def _compute_core_num(G):
        n = G.vertex_num
        degrees = [len(G.adj_list[u]) for u in range(n)]
        max_deg = max(degrees) if degrees else 0
        bucket = [[] for _ in range(max_deg+1)]
        for u in range(n):
            bucket[degrees[u]].append(u)

        core_num = [0] * n
        curr_degrees = degrees[:]

        for d in range(max_deg+1):
            while bucket[d]:
                u = bucket[d].pop()
                core_num[u] = d
                for w in G.adj_list[u]:
                    if curr_degrees[w] > d:
                        old = curr_degrees[w]
                        curr_degrees[w] -= 1
                        bucket[old].remove(w)
                        bucket[old-1].append(w)
        return core_num

    @staticmethod
    def _compute_nei_subgraph(G, neighbors_v, k):
      adj_sub = {}
      deg_sub = {}
      neighbors_set = set(neighbors_v)
      for u in neighbors_v:
          neighbours_in_sub = [w for w in G.adj_list[u] if w in neighbors_set]
          adj_sub[u] = neighbours_in_sub
          deg_sub[u] = len(neighbours_in_sub)

      peel_queue = deque(u for u in neighbors_v if deg_sub[u] < k)
      removed_set = set()
      while peel_queue:
          u0 = peel_queue.popleft()
          if u0 in removed_set:
              continue
          removed_set.add(u0)
          for w in adj_sub[u0]:
              if w not in removed_set:
                  deg_sub[w] -= 1
                  if deg_sub[w] < k:
                      peel_queue.append(w)

      core_vertices = [u for u in neighbors_v if u not in removed_set]

      tau_value = 0
      visited = set()
      for u in core_vertices:
          if u in visited:
              continue
          tau_value += 1
          stack = [u]
          visited.add(u)
          while stack:
              x = stack.pop()
              for w in adj_sub[x]:

                  if w not in removed_set and w not in visited:
                      visited.add(w)
                      stack.append(w)

      return tau_value

    @staticmethod
    def process(G, k):

        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # sd[v] -> τ_k(v) for all v
        """
        core_num = kCoreBaseStructuralDiversity._compute_core_num(G)

        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors_v = [u for u in G.adj_list[v] if core_num[u] >= k+1]
            if not neighbors_v:
                sd[v] = 0
            else:
                sd[v] = kCoreBaseStructuralDiversity._compute_nei_subgraph(G, neighbors_v, k)

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
