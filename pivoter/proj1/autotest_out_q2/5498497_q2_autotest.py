#!/usr/bin/env python3
# Auto-generated for 5498497

STUDENT_ID = "5498497"
STUDENT_NAME = "Kaishu Zhang"

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
        n = G.vertex_num #num of ertex
        sd = [0] * n #initialize result list

        for v in range(n):

          neighbors = G.adj_list[v] #get v neighbor list
          if not neighbors:
            continue

          #create sub_G - induced subgraph of the neighbor set
          sub_G = kCoreBaseStructuralDiversity._build_subgraph(G, neighbors)

          #compute num of k-core connected components
          sd[v] = kCoreBaseStructuralDiversity._compute_k_core(sub_G, k)

        return sd

    #Construct the induced subgraph
    @staticmethod
    def _build_subgraph(G, neighbors):
      sub_G = {u: set() for u in neighbors}
      for u in neighbors:
        for w in G.adj_list[u]:
          if w in neighbors:
            sub_G[u].add(w)

      return sub_G

    #k-core decomposition and count the connected components
    @staticmethod
    def _compute_k_core(sub_G, k):
      deg = {u: len(sub_G[u]) for u in sub_G}
      deleted = set()

      initial_nodes = []
      for u in sub_G:
        if deg[u] < k:
          initial_nodes.append(u)
      queue = deque(initial_nodes)

      while queue:
        u = queue.popleft()
        if u in deleted:
          continue

        deleted.add(u)
        for v in sub_G[u]:
          if v not in deleted:
            deg[v] -= 1
            if deg[v] < k:
              queue.append(v)

      remains = []
      for u in sub_G:
        if u not in deleted:
          remains.append(u)

      if not remains:
        return 0

      visited = set()
      cnt = 0
      for node in remains:
        if node not in visited:
          kCoreBaseStructuralDiversity._dfs_kcore(sub_G, node, visited, deleted)
          cnt += 1

      return cnt

    #dfs on the remaining nodes to count the num of connected components
    @staticmethod
    def _dfs_kcore(graph, start, visited, deleted):
      stack = [start]
      visited.add(start)

      while stack:
        u = stack.pop()
        for v in graph[u]:
          if v not in visited and v not in deleted:
            visited.add(v)
            stack.append(v)


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
