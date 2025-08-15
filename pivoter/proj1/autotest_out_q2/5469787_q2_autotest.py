#!/usr/bin/env python3
# Auto-generated for 5469787

STUDENT_ID = "5469787"
STUDENT_NAME = "Yiming Hou"

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

        def k_core_decomposition(sub_adj, k):
            """
            This function receives the adjacency list and k, used to filter out vertices with a deg less than k, and repeatedly excludes them,
            finally returns a set of vertices with degree greater than or equal to k.
            """
            # Filter out the list q where deg is less than k
            deg = {u: len(sub_adj[u]) for u in sub_adj}
            q = deque([u for u in sub_adj if deg[u] < k])
            # record what points were deleted in total to avoid duplication
            in_queue = set(q)
            while q:
                u = q.popleft()
                for v in sub_adj[u]:
                    # Reduce the deg of the point v adjacent to the deleted point u, and determine whether v should also be deleted after the reduction
                    if v in deg:
                        deg[v] -= 1
                        if deg[v] < k and v not in in_queue:
                            q.append(v)
                            in_queue.add(v)
                # Delete the point u whose deg is less than k
                if u in deg:
                    del deg[u]
            return set(deg.keys())  # Remaining nodes in k-core

        def get_connected_components(sub_adj, nodes):
            """
            This function receives the adjacency list and vertex list, and returns the connected graph composed of these vertices.
            """
            visited = set()     # Used to record visited nodes to prevent repeated traversal
            components = []     # Used to store all connected components (each component is a set)
            for node in nodes:
                if node not in visited:
                    comp = []       # The node set of the current connected component
                    stack = [node]     # The stack used by DFS
                    while stack:
                        u = stack.pop()
                        if u not in visited:
                            visited.add(u)
                            comp.append(u)
                            # Traverse adjacent nodes and only expand neighbors that are in the nodes set and have not been visited
                            for v in sub_adj[u]:
                                if v in nodes and v not in visited:
                                    stack.append(v)
                    components.append(set(comp))
            return components

        # Traverse each node v and calculate its structural diversity value sd[v]
        for v in range(n):
          neighbors = G.adj_list[v]
          # Edge cases: isolated vertices, directly set sd[v] = 0 and skip
          if neighbors==[]:
              sd[v] = 0
              continue
          # Construct a subgraph consisting of v's neighbors (excluding v itself)
          # Edge cases: graphs with no edges If there is no edge, neighbors will become (), sub_adj will become {},
          # then final_adj is also empty, and sd[v] will be equal to 0
          neighbors = set(neighbors)
          sub_adj = {u: [] for u in neighbors}

          # Build the adjacency list of the subgraph, keeping only the edges between neighbors
          for u in neighbors:
              for w in G.adj_list[u]:
                  if w in neighbors:
                      sub_adj[u].append(w)
          # use k_core_decomposition to exclude the points with insufficient degree
          k_core_nodes = k_core_decomposition(sub_adj, k)

          # construct the subgraph final_adj composed of vertices with degree greater than or equal to k
          final_adj = {u: [] for u in k_core_nodes}
          for u in k_core_nodes:
              for w in G.adj_list[u]:
                  if w in k_core_nodes:
                      final_adj[u].append(w)
          # Use get_connected_components to calculate the connected graph list, and directly assign the number of len to sd [v],
          components = get_connected_components(final_adj, k_core_nodes)
          sd[v] = len(components)
        # Return the structural diversity results of all nodes
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
