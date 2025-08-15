#!/usr/bin/env python3
# Auto-generated for 5505856

STUDENT_ID = "5505856"
STUDENT_NAME = "Grace Gao"

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
        sd = [0] * n

        for v in range(n):
          neighbors = G.adj_list[v]

          # Early termination for impossible cases
          if len(neighbors) < k:
              sd[v] = 0
              continue

          # Build neighbor-induced subgraph
          neighbor_adj = kCoreBaseStructuralDiversity._build_neighbor_subgraph(G, v)

          # Find k-cores and count components
          if k == 0:
              # For k=0, all vertices are in 0-cores
              sd[v] = kCoreBaseStructuralDiversity._count_connected_components(neighbor_adj)
          else:
              k_cores = kCoreBaseStructuralDiversity._find_k_cores(neighbor_adj, k)
              sd[v] = kCoreBaseStructuralDiversity._count_connected_components(k_cores)

        return sd

    @staticmethod
    def _build_neighbor_subgraph(G, v):
        """Build neighbor-induced subgraph for vertex v."""
        neighbors = G.adj_list[v]
        if not neighbors:
            return []

        neighbor_set = set(neighbors)
        neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbors)}

        neighbor_adj = [[] for _ in range(len(neighbors))]
        for i, u in enumerate(neighbors):
            for w in G.adj_list[u]:
                if w in neighbor_set and w != u:
                    neighbor_adj[i].append(neighbor_to_idx[w])

        return neighbor_adj

    @staticmethod
    def _find_k_cores(adj_list, k):
        """
        Find k-core subgraph using degeneracy ordering algorithm.
        Iteratively removes vertices with degree < k.
        """
        n = len(adj_list)
        if n == 0:
            return []

        degrees = [len(adj_list[i]) for i in range(n)]
        removed = [False] * n
        queue = deque()

        # Initialize queue with vertices having degree < k
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)

        # Iteratively remove vertices with degree < k
        while queue:
            u = queue.popleft()
            if removed[u]:
                continue

            removed[u] = True

            # Update neighbors' degrees
            for v in adj_list[u]:
                if not removed[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)

        # Build k-core subgraph
        vertex_mapping = {}
        new_id = 0

        for i in range(n):
            if not removed[i]:
                vertex_mapping[i] = new_id
                new_id += 1

        k_core_adj = []
        for i in range(n):
            if not removed[i]:
                new_adj = []
                for j in adj_list[i]:
                    if not removed[j]:
                        new_adj.append(vertex_mapping[j])
                k_core_adj.append(new_adj)

        return k_core_adj

    @staticmethod
    def _count_connected_components(adj_list):
        """Count connected components using iterative DFS."""
        n = len(adj_list)
        if n == 0:
            return 0

        visited = [False] * n
        components = 0

        for i in range(n):
            if not visited[i]:
                # Iterative DFS to avoid stack overflow
                stack = [i]
                visited[i] = True

                while stack:
                    u = stack.pop()
                    for v in adj_list[u]:
                        if not visited[v]:
                            visited[v] = True
                            stack.append(v)

                components += 1

        return components



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
