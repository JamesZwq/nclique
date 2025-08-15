#!/usr/bin/env python3
# Auto-generated for 5511407

STUDENT_ID = "5511407"
STUDENT_NAME = "Cheney Wang"

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
        #return sd


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
        for i in range(n):
          # get the adjacent list of vertex i
          nbrs = G.adj_list[i]
          if nbrs:
              nbrset = set(nbrs)
              # initial
              induced_subgraph = {}
              # search any adjacent vertex x of i
              for x in nbrs:
                  # initial
                  induced_subgraph[x] = []
                  # search any adjacent vertex y of x
                  for y in G.adj_list[x]:
                      if y in nbrset:
                          induced_subgraph[x].append(y)
              # calculate kcore
              sd[i] = kCoreBaseStructuralDiversity.Kcore(induced_subgraph, k)
          else:
              sd[i] = 0
        return sd

    @staticmethod
    def Kcore(newg, k):
      if not newg:
          return 0
      #calculate degree of vertex
      deg = {i: len(adj) for i, adj in newg.items()}
      #remove the degree that less than k
      removed = set(i for i, d in deg.items() if d < k)
      #using double end queue restore vertex,which easy to out and end
      rq = deque(removed)
      concomp = 0
      stack = []
      #loop to delete the vertex that less than k
      while rq:
        i = rq.popleft()
        #remove the adjacent vertex of i
        for j in newg[i]:
          if j not in removed:
            #degree minus 1
              deg[j] -= 1
              if deg[j] < k:
                  removed.add(j)
                  rq.append(j)
      #the vertex not removed
      P = set(i for i in newg if i not in removed)
      if not P:
          return 0
      #calculate connected components
      while P:
          #connected component add 1
          concomp += 1
          #select 1 as the start
          starpoint = P.pop()
          stack.append(starpoint)
          #using DFS,remove vertex
          while stack:
              i = stack.pop()
              for j in newg[i]:
                  if j in P:
                      P.remove(j)
                      stack.append(j)
      return concomp

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
