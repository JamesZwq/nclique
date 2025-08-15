#!/usr/bin/env python3
# Auto-generated for 5521872

STUDENT_ID = "5521872"
STUDENT_NAME = "Ruochen Zhang"

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
        adj = G.adj_list
        core = kCoreBaseStructuralDiversity.CoreDecomposition(adj)

        result = [0] * n
        valid = [0] * n
        alive = [0] * n
        visited = [0] * n
        dloc = [0] * n

        for v in range(n):
          verifier = v + 1
          neighbor = []
          subgraph_cnt = 0
          for u in adj[v]:
            if core[u]:
                valid[u] = verifier
                neighbor.append(u)

          if not neighbor:
            result[v] = 0
            continue

          for u in neighbor:
            cnt = 0
            for w in adj[u]:
              if valid[w] == verifier:
                cnt += 1
            dloc[u] = cnt
            alive[u] = verifier

          if k > 0:
            out_buffer = deque([u for u in neighbor if dloc[u] < k])
            while out_buffer:
              u = out_buffer.popleft()
              if alive[u] != verifier:
                continue
              alive[u] = 0
              for w in adj[u]:
                if valid[w] == verifier and alive[w] == verifier:
                  dloc[w] -= 1
                  if dloc[w] == k-1:
                    out_buffer.append(w)
          for u in neighbor:
            if alive[u] == verifier and visited[u] != verifier:
              visited[u] = verifier
              subgraph_cnt += 1
              q = deque([u])
              while q:
                temp = q.popleft()
                for w in adj[temp]:
                  if valid[w] == verifier and alive[w] == verifier and visited[w] != verifier:
                    visited[w] = verifier
                    q.append(w)
          result[v] = subgraph_cnt
        return result


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def CoreDecomposition(adj):
      n = len(adj)
      if n == 0:
        return []

      deg = [len(adj[v]) for v in range(n)]
      maxd = max(deg) if n > 0 else 0
      bin_cnt = [0] * (maxd + 1)
      for d in deg:
        bin_cnt[d] += 1
      start = 0
      bin_pos = [0] * (maxd + 1)
      for d in range(maxd + 1):
        bin_pos[d] = start
        start += bin_cnt[d]
      vert = [0] * n
      pos = [0] * n
      next_pos = bin_pos[:]
      for v in range(n):
        d = deg[v]
        vert[next_pos[d]] = v
        pos[v] = next_pos[d]
        next_pos[d] += 1

      for i in range(n):
        v = vert[i]
        for u in adj[v]:
          if deg[u] > deg[v]:
            du = deg[u]
            pu = pos[u]
            pw = bin_pos[du]
            w = vert[pw]
            if u != w:
              vert[pu], vert[pw] = vert[pw], vert[pu]
              pos[u], pos[w] = pw, pu
            bin_pos[du] += 1
            deg[u] -= 1
      return deg

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
