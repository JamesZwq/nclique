#!/usr/bin/env python3
# Auto-generated for 5452066

STUDENT_ID = "5452066"
STUDENT_NAME = "Rin Ohsugi"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def compute_core_numbers(adj):
        n = len(adj)
        deg = [len(adj[v]) for v in range(n)]
        max_deg = max(deg) if n > 0 else 0

        bin = [0] * (max_deg + 1)
        for d in deg:
          bin[d] += 1

        start = 0
        for d in range(max_deg + 1):
          count = bin[d]
          bin[d] = start
          start += count

        pos = [0] * n
        vert = [0] * n

        for v in range(n):
          d = deg[v]
          idx = bin[d]
          pos[v] = idx
          vert[idx] = v
          bin[d] += 1

        for d in range(max_deg, 0, -1):
          bin[d] = bin[d - 1]
        bin[0] = 0

        for i in range(n):
          v = vert[i]
          for u in adj[v]:
            if deg[u] > deg[v]:
              du = deg[u]
              pu = pos[u]
              pw = bin[du]
              w = vert[pw]

              if u != w:
                vert[pu], vert[pw] = vert[pw], vert[pu]
                pos[u], pos[w] = pw, pu

              bin[du] += 1
              deg[u] -= 1

        return deg

    @staticmethod
    def get_connected_components(adj_list, nodes):
      visited = set()
      comps = []
      node_set = set(nodes)

      def dfs(u, comp):
        visited.add(u)
        comp.append(u)
        if len(adj_list) > u:
          for v in adj_list[u]:
            if v in node_set and v not in visited:
              dfs(v, comp)

      for u in nodes:
        if u not in visited:
          comp = []
          dfs(u, comp)
          comps.append(comp)

      return comps


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
        result = [0] * n

        for v in range(n):
            nbrs = G.adj_list[v]
            if not nbrs:
                continue

            index_map = {u: i for i, u in enumerate(nbrs)}
            sub_adj = [[] for _ in range(len(nbrs))]

            for i, u in enumerate(nbrs):
                for w in G.adj_list[u]:
                    if w in index_map:
                        sub_adj[i].append(index_map[w])

            core_nums = kCoreBaseStructuralDiversity.compute_core_numbers(sub_adj)

            valid_nodes = [i for i, c in enumerate(core_nums) if c >= k]
            components = kCoreBaseStructuralDiversity.get_connected_components(sub_adj, valid_nodes)

            result[v] = len(components)

        return result


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
