#!/usr/bin/env python3
# Auto-generated for 5602082

STUDENT_ID = "5602082"
STUDENT_NAME = "Yizhu Zhou"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        result = []

        # for every vertiex in the graph
        for v in range(n):
            # get its neighbour set
            nbr_nodes = set(G.adj_list[v])

            # get its neighbour-induced subgraph by inputing its neighbour set
            subgraph = kCoreBaseStructuralDiversity.get_induced_subgraph(G, nbr_nodes)

            # Repeatedly delete the nodes with degrees less than k to get the k-core nodes of the graph
            k_core_nodes = kCoreBaseStructuralDiversity.k_core_decomposition(subgraph, k)

            num_components = kCoreBaseStructuralDiversity.count_components(subgraph, k_core_nodes)

            result.append(num_components)

        return result

    @staticmethod
    def get_induced_subgraph(G, nodes):
        subgraph = {}
        for u in nodes:
            subgraph[u] = []    # 2: [0, 1, 3]

        # for each node u
        for u in nodes:
            # travel through all veritces in the graph to sort out their edges
            for v in G.adj_list[u]:
                if v in nodes:
                    subgraph[u].append(v)

        return subgraph



    @staticmethod
    def k_core_decomposition(subgraph, k):

        # count the degree of each nodes
        degrees = {}
        for node in subgraph:
            degrees[node] = len(subgraph[node])

        # record all the nodes with degree < k
        queue = deque()
        for node in degrees:
            if degrees[node] < k:
                queue.append(node)

        removed_nodes = set()

        # remove nodes from the queue and update the degrees of neighbors
        while queue:
            current = queue.popleft()
            removed_nodes.add(current)

            for neighbor in subgraph[current]:

                # skip the node if it has been removed
                if neighbor in removed_nodes:
                    continue
                degrees[neighbor] -= 1

                if degrees[neighbor] == k - 1:
                    queue.append(neighbor)

        # return all the remaining nodes(k-core nodes)
        initial_nodes = set(subgraph.keys())
        return initial_nodes - removed_nodes

    @staticmethod   # use BFS to count the number of connected component using the result of the remaining k-core nodes
    def count_components(subgraph, valid_nodes):
        visited = set()
        count = 0

        for node in valid_nodes:
          if node in visited:
              continue

          count += 1
          queue = deque()
          queue.append(node)
          visited.add(node)

          while queue:
              current = queue.popleft()
              # travel through the neighbours of current node
              for neighbor in subgraph[current]:
                  if neighbor in valid_nodes and neighbor not in visited:
                      visited.add(neighbor)
                      queue.append(neighbor)

        return count

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
