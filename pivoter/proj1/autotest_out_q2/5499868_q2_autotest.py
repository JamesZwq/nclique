#!/usr/bin/env python3
# Auto-generated for 5499868

STUDENT_ID = "5499868"
STUDENT_NAME = "Chufan Wu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # Build the neighbor-induced subgraph for a given vertex.
    @staticmethod
    def build_neighbor_induced_subgraph(G, v):
      neighbors = G.adj_list[v]
      nbr_set = set(neighbors)
      sub_adj = defaultdict(list)
      for u in neighbors:
          for w in G.adj_list[u]:
              if w in nbr_set:
                  sub_adj[u].append(w)
      return sub_adj, nbr_set

    # Count connected components among a subset of nodes(alive set) in a subgraph.
    @staticmethod
    def count_components_alive(sub_adj, alive):
        visited = set()
        comps = 0
        for s in alive:
            if s in visited:
                continue
            comps += 1
            q = deque([s])
            visited.add(s)
            while q:
                u = q.popleft()
                for w in sub_adj[u]:
                    if w in alive and w not in visited:
                        visited.add(w)
                        q.append(w)
        return comps

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            # No neighbors means zero components
            if not neighbors:
                sd[v] = 0
                continue

            # Build neighbor-induced subgraph
            sub_adj, nbr_set = kCoreBaseStructuralDiversity.build_neighbor_induced_subgraph(G,v)

            # Initialize alive set and degree map
            alive = set(nbr_set)
            deg = {u: 0 for u in nbr_set}
            for u in nbr_set:
                deg[u] = sum(1 for w in sub_adj[u] if w in nbr_set)

            # if every neighbor has degree >= k, count components directly
            if all(d >= k for d in deg.values()):
                sd[v] = kCoreBaseStructuralDiversity.count_components_alive(sub_adj, alive)
                continue

            # k-core peeling: iteratively remove nodes with degree < k
            q = deque([u for u in nbr_set if deg[u] < k])
            while q:
                x = q.popleft()
                if x not in alive:
                    continue
                alive.remove(x)
                for y in sub_adj[x]:
                    if y in alive:
                        deg[y] -= 1
                        if deg[y] < k:
                            q.append(y)

            # Count connected components among remaining alive nodes
            sd[v] = kCoreBaseStructuralDiversity.count_components_alive(sub_adj, alive)

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
