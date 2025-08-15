#!/usr/bin/env python3
# Auto-generated for 5495976

STUDENT_ID = "5495976"
STUDENT_NAME = "Hanwen Miao"

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

        degree = []
        for v in range(n):
            deg = len(G.adj_list[v])
            degree.append(deg)

        if n == 0:
            return sd

        if k < 0:
            return sd

        core_level, core_k1_nodes = kCoreBaseStructuralDiversity._find_max_k_core_per_node(G, degree, k)
        if not core_k1_nodes:
            return sd

        components = kCoreBaseStructuralDiversity._find_components_bfs(G, core_k1_nodes)
        for component in components:
          connected_v = set()
          for v in component:
              connected_v.add(v)
              connected_v.update(G.adj_list[v])

          for v in connected_v:
              kCoreBaseStructuralDiversity._process_vertex_tau_k(
                  G, v, k, core_level, sd
              )

        return sd


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def _sort_vertices_by_degree(degree):
        if not degree:
            return []

        max_deg = max(degree)
        buckets = [[] for _ in range(max_deg + 1)]

        for v, d in enumerate(degree):
            buckets[d].append(v)

        sorted_vertices = []
        for bucket in buckets:
            sorted_vertices.extend(bucket)

        return sorted_vertices

    @staticmethod
    def _find_max_k_core_per_node(G, degree, k):
        n = len(degree)
        if n == 0:
            return [], set()

        core_level = [0] * n
        sorted_vertices = kCoreBaseStructuralDiversity._sort_vertices_by_degree(degree)
        deg = degree.copy()

        for i in range(n):
            v = sorted_vertices[i]
            core_level[v] = deg[v]
            for u in G.adj_list[v]:
                if deg[u] > deg[v]:
                    deg[u] -= 1

        core_k1_nodes = set()
        for v in range(n):
            if core_level[v] >= k + 1:
                core_k1_nodes.add(v)

        return core_level, core_k1_nodes

    @staticmethod
    def _find_components_bfs(G, node_subset):
      visited_set = set()
      result = []

      for start_node in node_subset:
          if start_node in visited_set:
              continue

          group = []
          queue = deque([start_node])
          visited_set.add(start_node)

          while queue:
              current = queue.popleft()
              group.append(current)

              for neighbor in G.adj_list[current]:
                  if neighbor in node_subset and neighbor not in visited_set:
                      visited_set.add(neighbor)
                      queue.append(neighbor)

          result.append(group)

      return result



    @staticmethod
    def _process_vertex_tau_k(G, v, k, core_level, sd):
        k_core_neighbors = []
        degrees_in_subgraph = {}

        for u in G.adj_list[v]:
            if core_level[u] >= k:
                k_core_neighbors.append(u)

        if len(k_core_neighbors) < k:
            return


        neighbors_set = set(k_core_neighbors)
        for u in k_core_neighbors:
            degree = sum(1 for w in G.adj_list[u] if w in neighbors_set)
            degrees_in_subgraph[u] = degree


        tau_k = kCoreBaseStructuralDiversity._compute_k_core_count(
            G, neighbors_set, k, degrees_in_subgraph
        )

        if tau_k > sd[v]:
          sd[v] = tau_k



    @staticmethod
    def _compute_k_core_count(G, vertex_set, k, precomputed_degrees):
        remain_nodes = kCoreBaseStructuralDiversity._peel_k_core(G, vertex_set, k, precomputed_degrees)

        if len(remain_nodes) == 0:
            return 0

        components = kCoreBaseStructuralDiversity._find_components_bfs(G, remain_nodes)
        return len(components)


    @staticmethod
    def _peel_k_core(graph, candidate_nodes, k, degree_map):
        current_degree = degree_map.copy()
        to_remove = deque()
        removed_nodes = set()

        for node in candidate_nodes:
            if current_degree[node] < k:
                to_remove.append(node)
                removed_nodes.add(node)

        while to_remove:
            node = to_remove.popleft()
            for neighbor in graph.adj_list[node]:
                if neighbor in candidate_nodes and neighbor not in removed_nodes:
                    current_degree[neighbor] -= 1
                    if current_degree[neighbor] < k:
                        to_remove.append(neighbor)
                        removed_nodes.add(neighbor)

        result = set()
        for node in candidate_nodes:
            if node not in removed_nodes:
                result.add(node)

        return result

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
