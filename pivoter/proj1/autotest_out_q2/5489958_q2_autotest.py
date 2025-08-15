#!/usr/bin/env python3
# Auto-generated for 5489958

STUDENT_ID = "5489958"
STUDENT_NAME = "Yifan Teng"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass
    # build subgraph for one vertex
    @staticmethod
    def build_subgraph(G,neighbors):
      subgraph={}
      for neighbor in neighbors:
        new_neighbors=set()
        for element in G.adj_list[neighbor]:
          if element in neighbors and element>neighbor:
            new_neighbors.add(element)
            if element in subgraph:
              subgraph[element].add(neighbor)
            else:
              subgraph[element]={neighbor}
        if neighbor not in subgraph:
          subgraph[neighbor]=new_neighbors
        else:
          subgraph[neighbor].update(new_neighbors)
      return subgraph


    # check if any vertex in the subgraph has a degree less than k.
    # return the vertex a vertex with degree less than k; return -1 if none.
    @staticmethod
    def check_subgraph(subgraph,k):
      for key,value in subgraph.items():
        if len(value)<k:
          return key
      return -1

    # delete every vertex whose degree is less than k
    @staticmethod
    def delete_vertices(graph,k):
      subgraph=graph
      delete_vertex=kCoreBaseStructuralDiversity.check_subgraph(subgraph,k)
      while delete_vertex!=-1:
        del subgraph[delete_vertex]
        for _,value in subgraph.items():
          if delete_vertex in value:
            value.remove(delete_vertex)
        delete_vertex=kCoreBaseStructuralDiversity.check_subgraph(subgraph,k)
      return subgraph

    # count the number of k-cores in the neighbour-induced subgraph
    # this function is based on the demo code in topic 1.1, with slight modifications
    @staticmethod
    def compute_kcore_number(subgraph):
      count=0
      visited={vertex: False for vertex in subgraph}
      for key,value in subgraph.items():
        if visited[key]==True:
          continue
        else:
          count+=1
          visited[key]=True
        if not value:
          continue
        queue=deque()
        queue.append(key)
        while len(queue)>0:
          s=queue.popleft()
          for neighbor in subgraph[s]:
            if visited[neighbor]==True:
              continue
            else:
              queue.append(neighbor)
              visited[neighbor]=True
      return count

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
        for i in range(0,n):
          neighbors=set(G.adj_list[i])
          if len(neighbors)==0:
            continue
          # build subgraph for one vertex
          subgraph=kCoreBaseStructuralDiversity.build_subgraph(G,neighbors)
          # delete every vertex whose degree is less than k
          final_subgraph=kCoreBaseStructuralDiversity.delete_vertices(subgraph,k)
          if not final_subgraph:
            continue
          # count the number of k-core
          else:
            count=kCoreBaseStructuralDiversity.compute_kcore_number(final_subgraph)
            sd[i]=count
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
