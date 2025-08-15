#!/usr/bin/env python3
# Auto-generated for 5527498

STUDENT_ID = "5527498"
STUDENT_NAME = "Zihan Chen"

# ======= 学生代码 =======
from inspect import stack
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
          nbr_v = kCoreBaseStructuralDiversity.get_nbr(G, v)      # induced subgraph with neighbors of v
          core_nodes = kCoreBaseStructuralDiversity.k_core(nbr_v, k)      # set of nodes satisfy the k-core condition
          cc = kCoreBaseStructuralDiversity.get_cc(core_nodes, nbr_v)   # get connected components
          sd[v] = len(cc)     # number of components
        return sd

    @staticmethod
    def get_nbr(G, v):
      nbr = G.adj_list[v]     # get neighbor of v
      nbr_set = set(nbr)     # deduplication
      nbr_adj = {a: [] for a in nbr_set}  # initialize the adjacency list
      for a in nbr_set:
        for b in G.adj_list[a]:     # b is a's neighbor
          if b in nbr_set:
            nbr_adj[a].append(b)    # add in adj list
      return nbr_adj

    @staticmethod
    def k_core(adj_list, k):
      deg = {a: len(adj_list[a]) for a in adj_list}  # initial degree list
      queue = deque()
      for a in deg:       # nodes that degree < k add to list
        if deg[a] < k:
          queue.append(a)
      while queue:
        a = queue.popleft()     # remove a node not sufficiently connected
        for b in adj_list[a]:
          if b in deg:
            deg[b] -= 1     # v degree-1
            if deg[b] == k - 1:
              queue.append(b)    # v not sufficiently connected, add to queue
        del deg[a]     # delete u
      return set(deg.keys())    # remaining nodes to be k-core

    @staticmethod
    def get_cc(nodes, adj_list):
      visited = set()
      cc = []

      for node in nodes:
        if node not in visited:
          current_cc = []     # current connected component
          stack = [node]
          while stack:        # DFS
            a = stack.pop()
            if a not in visited:
              visited.add(a)
              current_cc.append(a)
              for b in adj_list[a]:
                if b not in visited and b in nodes:
                  stack.append(b)
          cc.append(current_cc)   # save current_cc
      return cc

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
