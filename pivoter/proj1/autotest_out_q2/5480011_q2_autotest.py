#!/usr/bin/env python3
# Auto-generated for 5480011

STUDENT_ID = "5480011"
STUDENT_NAME = "Justin Xia"

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

        def dfs(u):
          if visited[u]:
            return

          visited[u] = True
          for w in G.adj_list[u]:
            if is_neighbour[w] and degree[w] >= k:
              dfs(w)

        n = G.vertex_num
        sd = [0] * n

        # Loop runs n times
        for v in range(n):

          # O(n) initialisation
          is_neighbour = [False] * n
          for u in G.adj_list[v]:
            is_neighbour[u] = True

          degree = [0] * n
          queue = deque()
          # O(n + m) to iterate through all vertices and edges
          for u in range(n):
            degree[u] = len([w for w in G.adj_list[u] if is_neighbour[w]])
            if degree[u] < k and is_neighbour[u]:
              queue.append(u)

          # Modified multi-source bfs, taking O(n + m)
          visited = [False] * n
          while queue:
            u = queue.popleft()
            if visited[u]:
              continue
            visited[u] = True

            for w in G.adj_list[u]:
              if not is_neighbour[w]:
                continue

              degree[w] -= 1
              if degree[w] < k:
                queue.append(w)

          # O(n) loop
          kcore = []
          for u in range(n):
            if degree[u] >= k and is_neighbour[u]:
              kcore.append(u)

          # Since we use a global visited array for every dfs, even though
          # we call dfs multiple times it will still only have time complexity
          # of O(n + m)
          visited = [False] * n
          for u in kcore:
            if visited[u]:
              continue
            dfs(u)
            sd[v] += 1

        # Work inside the major loop does O(n + m) work, and the loop runs n times,
        # so the final time complexity is O(n(n + m))

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
