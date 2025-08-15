#!/usr/bin/env python3
# Auto-generated for 5514271

STUDENT_ID = "5514271"
STUDENT_NAME = "Han Dong"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    def find_neighor_induced_subgraph(G,current_node):
      neighbors = set(G.adj_list[current_node])
      #use dict to store neighbor induceed subgraph
      pre_subgraph = {}
      for neighbor in neighbors:         #O(|V|)
        new_neighbors = []
        for node in G.adj_list[neighbor]:   #O(|V|)
          if node in neighbors:
            new_neighbors.append(node)
        pre_subgraph[neighbor] = new_neighbors
      return pre_subgraph

    def DFS(current_subgraph,start,visited):
      visited.add(start)
      for neighbor in current_subgraph[start]:
        if neighbor not in visited:
          kCoreBaseStructuralDiversity.DFS(current_subgraph,neighbor,visited)

    def comuptue_num_connected_components(current_subgraph):
      #use DFS to compute the number of connected components
      #record visited node
      visited = set()
      num = 0
      for node in current_subgraph.keys():
        if node not in visited:
          kCoreBaseStructuralDiversity.DFS(current_subgraph,node,visited)
          num += 1
      return num

    def comupte_num_k_core(current_subgraph, k):
      queue = deque()
      degree = {}
      removed = set()
      #find nodes whose degree less than k and push it in queue and intialize a dict to store the degrees of each node
      for node in current_subgraph.keys():        #O(|v|)
        degree[node] = len(current_subgraph[node])
        if degree[node] < k:
          queue.append(node)
          removed.add(node)
      #Iteratively remove the nodes whose degree less than k
      #O(|v|+|e|)
      while queue:
        del_node = queue.popleft()
        for neighbor in current_subgraph[del_node]:
          if neighbor not in removed:
            degree[neighbor] -= 1
            if degree[neighbor] < k:
              queue.append(neighbor)
              removed.add(neighbor)
      #create new subgraph contain the rest of nodes
      #O(|v|+|e|)
      new_subgraph = {}
      for node in current_subgraph.keys():
        if node not in removed:
          new_neighbors = []
          for neighbor in current_subgraph[node]:
            if neighbor not in removed:
              new_neighbors.append(neighbor)
          new_subgraph[node] = new_neighbors
      #use DFS to compute the number of connected components
      #O（|v|+|e|）
      return kCoreBaseStructuralDiversity.comuptue_num_connected_components(new_subgraph)

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

        for current_node in range(n):
          #get neighbor induced subgraph
          current_subgraph = kCoreBaseStructuralDiversity.find_neighor_induced_subgraph(G,current_node)  #O(|V|*|V|)
          #compute k cores
          sd[current_node] = kCoreBaseStructuralDiversity.comupte_num_k_core(current_subgraph,k)  #O（|V|*|V|）

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
