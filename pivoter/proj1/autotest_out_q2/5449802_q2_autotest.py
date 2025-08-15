#!/usr/bin/env python3
# Auto-generated for 5449802

STUDENT_ID = "5449802"
STUDENT_NAME = "Wencheng Guan"

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
        # all the vertecies start from 0
        for i in range(n):
          neighbour=set(G.adj_list[i])
          # make sure this spot is not neighbour and haven't been visited过
          visited=set()
          group=set()
          # we use DFS to find the neighbour
          stack=deque()
          component=[]
          for l in neighbour:
            if l not in visited:
              visited.add(l)
              stack.append(l)
              group=set([l])
            else:
              continue
            while stack:
              v=stack.pop()
              for j in G.adj_list[v]:
                if j in neighbour and j not in visited:
                  visited.add(j)
                  stack.append(j)
                  group.add(j)
            component.append(group)
          # from here, we already find all the subgraph
          number=0
          for subset in component:
            number+=kCoreBaseStructuralDiversity.k_core(G,subset,k)
          sd[i]=number

        return sd
    def k_core(G,subgraph,k):


      clean_graph = {
          u: [v for v in G.adj_list[u] if v in subgraph]
          for u in subgraph
      }
      # using dictionary to store each vertex's degree
      degree = {u: len(neigh) for u, neigh in clean_graph.items()}
      delete=set()
      queue = deque()
      for u in clean_graph:
        if degree[u]<k:
          queue.append(u)
          delete.add(u)

      # while queue can make sure we will keep delete the vertex until the queue is empty

      while queue:
        # get the node from the bottom
        node=queue.popleft()
        delete.add(node)
        for nei in clean_graph[node]:
          if nei in degree and nei not in delete:
            degree[nei]-=1
            if degree[nei]<k:
              queue.append(nei)
              delete.add(nei)
        # make sure we delete this vertex

        del degree[node]
      # these are the vertices haven't been deleted
      remain=set(degree.keys())
      if len(remain)==0:
        return 0
      # after this, see h
      visited=set()
      number=0
      line=deque()
      for i in remain:
        if i not in visited:

          number+=1
          visited.add(i)
          line.append(i)
          while line:
            v=line.popleft()
            for j in G.adj_list[v]:
              if j in remain and j not in visited and j not in delete:
                visited.add(j)
                line.append(j)

      return number


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
