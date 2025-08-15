#!/usr/bin/env python3
# Auto-generated for 5522711

STUDENT_ID = "5522711"
STUDENT_NAME = "Haolin Huang"

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
        n = G.vertex_num
        result = []

        for v in range(n):
          #find all
          neighbors = []
          for neighbor in G.adj_list[v]:
              neighbors.append(neighbor)

          if len(neighbors) == 0:
              result.append(0)
              continue

          # buuld Neighbour-induced subgraph
          neighbor_graph = {}
          for n in neighbors:
              neighbor_graph[n] = []

          for n in neighbors:
              for v in G.adj_list[n]:
                  if v in neighbors:
                      neighbor_graph[n].append(v)

          # find all nodes belong k-core
          k_core_nodes = kCoreBaseStructuralDiversity.find_k_core_nodes(neighbor_graph, k)

          if len(k_core_nodes) == 0:
              result.append(0)
              continue

          #count k core
          count = kCoreBaseStructuralDiversity.count_k_core(neighbor_graph, k_core_nodes)
          result.append(count)

        return result

    @staticmethod
    def find_k_core_nodes(graph, k):
        if k == 0:   #if k=0 return all
            return list(graph.keys())


        temp_graph = {}  #copy graph
        degrees = {}

        for node in graph:
            temp_graph[node] = []
            for neighbor in graph[node]:
                temp_graph[node].append(neighbor)
            degrees[node] = len(temp_graph[node])


        nodes_to_remove = []
        for node in temp_graph:
            if degrees[node] < k:
                nodes_to_remove.append(node)

        removed_nodes = []
        #save node that in neighbor
        while len(nodes_to_remove) > 0:
            current_node = nodes_to_remove.pop(0)

            if current_node in removed_nodes:
                continue

            removed_nodes.append(current_node)

            for neighbor in temp_graph[current_node]:
                if neighbor not in removed_nodes:
                    degrees[neighbor] = degrees[neighbor] - 1
                    if degrees[neighbor] < k and neighbor not in nodes_to_remove:
                        nodes_to_remove.append(neighbor)

        k_core_nodes = []
        for n in temp_graph:
            if n not in removed_nodes:
                k_core_nodes.append(n)

        return k_core_nodes

    @staticmethod
    def count_k_core(graph, nodes):
        if len(nodes) == 0:
            return 0

        visited = []
        count = 0

        #dfs
        for n in nodes:
            if n not in visited:
                count = count + 1

                nodes_to_visit = [n]

                while len(nodes_to_visit) > 0:
                    current_node = nodes_to_visit.pop()

                    if current_node in visited:
                        continue

                    visited.append(current_node)

                    for neighbor in graph[current_node]:
                        if neighbor in nodes and neighbor not in visited:
                            nodes_to_visit.append(neighbor)

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
