#!/usr/bin/env python3
# Auto-generated for 5447603

STUDENT_ID = "5447603"
STUDENT_NAME = "Kailang Zhang"

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
       
       Time Complexity: O(n * m)
       where n is the number of vertices and m is the number of edges
       """
       n = G.vertex_num
       tau_k_values = [0] * n
       
       for v in range(n):
           neighbors = G.adj_list[v]
           
           # Optimization: if the number of neighbors is insufficient to form a k-core, skip
           # A non-empty k-core requires at least k+1 vertices
           if len(neighbors) <= k:
               tau_k_values[v] = 0
               continue
           
           # Convert neighbor list to set for efficient lookup
           neighbor_set = set(neighbors)
           
           # Build neighbor-induced subgraph
           # Map original vertex IDs to new IDs in subgraph (0, 1, 2, ...)
           vertex_map = {node_id: i for i, node_id in enumerate(neighbors)}
           subgraph_adj = [[] for _ in range(len(neighbors))]
           
           for u_original in neighbors:
               u_subgraph = vertex_map[u_original]
               for w_original in G.adj_list[u_original]:
                   # Check if w_original is also a neighbor of v
                   if w_original in neighbor_set:
                       w_subgraph = vertex_map[w_original]
                       subgraph_adj[u_subgraph].append(w_subgraph)
           
           # Count k-cores in the subgraph
           tau_k_values[v] = kCoreBaseStructuralDiversity._count_k_cores_in_subgraph(subgraph_adj, k)
           
       return tau_k_values


   ################################################################################
   # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
   ################################################################################
   
   @staticmethod
   def _count_k_cores_in_subgraph(subgraph_adj, k):
       """
       Efficiently count the number of k-cores in a subgraph.
       This function first performs k-core decomposition, then counts connected components in the remaining graph.
       
       Time Complexity: O(m_sub)
       where m_sub is the number of edges in the subgraph
       """
       n = len(subgraph_adj)
       if n == 0:
           return 0

       # 1. k-core decomposition
       degrees = [len(neighbors) for neighbors in subgraph_adj]
       q = deque()
       
       # is_removed[i] = True indicates vertex i is removed during k-core decomposition
       is_removed = [False] * n

       for i in range(n):
           if degrees[i] < k:
               q.append(i)
               is_removed[i] = True
       
       while q:
           u = q.popleft()
           for v_neighbor in subgraph_adj[u]:
               if not is_removed[v_neighbor]:
                   degrees[v_neighbor] -= 1
                   if degrees[v_neighbor] < k:
                       q.append(v_neighbor)
                       is_removed[v_neighbor] = True
       
       # 2. Count connected components in the remaining graph
       # After k-core decomposition, each connected component in the remaining subgraph is a k-core
       num_k_cores = 0
       visited = [False] * n
       
       for i in range(n):
           # If vertex i is not removed and not visited, it belongs to a new k-core
           if not is_removed[i] and not visited[i]:
               num_k_cores += 1
               component_q = deque([i])
               visited[i] = True
               
               while component_q:
                   u = component_q.popleft()
                   for v_neighbor in subgraph_adj[u]:
                       # Perform BFS among vertices that are also not removed
                       if not is_removed[v_neighbor] and not visited[v_neighbor]:
                           visited[v_neighbor] = True
                           component_q.append(v_neighbor)
       
       return num_k_cores

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
