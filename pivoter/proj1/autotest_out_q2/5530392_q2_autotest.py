#!/usr/bin/env python3
# Auto-generated for 5530392

STUDENT_ID = "5530392"
STUDENT_NAME = "Xiangting Chen"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # peel k‑core inside a small subgraph
    @staticmethod
    def _k_core_peel(sub_adj, k):
      m = len(sub_adj)
      degree = [len(lst) for lst in sub_adj]# degree of each vertex
      alive = [True] * m
      Q = deque(i for i, d in enumerate(degree) if d < k)# vertices to delete

      while Q:
        v = Q.popleft()
        if not alive[v]:
          continue
        alive[v] = False# remove v
        for u in sub_adj[v]:# update neighbours
          if alive[u]:
            degree[u] -= 1
            if degree[u] == k - 1:
              Q.append(u)
      return alive

    # main API
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
        tau = [0] * n# answer array

        # build neighbour induced subgraph
        def build_subgraph(neighbours):
          local_id = {u: idx for idx, u in enumerate(neighbours)}
          m = len(neighbours)
          sub_adj = [[] for _ in range(m)]
          for u in neighbours:
            for w in G.adj_list[u]:
              if w in local_id and w != u:
                sub_adj[local_id[u]].append(local_id[w])
          return sub_adj

        # main loop
        for v in range(n):
          nbrs = G.adj_list[v]# neighbours of v

          # low degree: τ=0
          if len(nbrs) < k:
            tau[v] = 0
            continue

          # neighbour induced subgraph
          sub_adj = build_subgraph(nbrs)

          # k‑core peeling in the subgraph
          alive = kCoreBaseStructuralDiversity._k_core_peel(sub_adj, k)

          # count CCs among alive vertices
          m = len(nbrs)
          seen = [False] * m
          comp_cnt = 0
          for i in range(m):
            if alive[i] and not seen[i]:
              comp_cnt += 1
              stack = [i]
              while stack:
                cur = stack.pop()
                if seen[cur]:
                  continue
                seen[cur] = True
                for nxt in sub_adj[cur]:
                  if alive[nxt] and not seen[nxt]:
                    stack.append(nxt)
          tau[v] = comp_cnt# τ_k(v)

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
