#!/usr/bin/env python3
# Auto-generated for 5485269

STUDENT_ID = "5485269"
STUDENT_NAME = "Pei-Yi Hsieh"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _peel_subgraph(num_nodes, adj_list, degrees, k):
        """
        Performs k-core decomposition using a queue-based peeling algorithm.
        This is a standard and efficient implementation.
        """
        q = deque()
        is_in_core = [True] * num_nodes

        # Initially, add all nodes with degree<k to removal queue
        for i in range(num_nodes):
            if degrees[i] < k:
                q.append(i)
                is_in_core[i] = False
        
        while q:
            u = q.popleft()
            for neighbor in adj_list[u]:
                # degree of a neighbor must be updated to reflect removal of u
                degrees[neighbor] -= 1
                
                # If the degree reduction causes neighbor to also fall below k,
                # and it has not been queued for removal before, add it to queue.
                if degrees[neighbor] < k and is_in_core[neighbor]:
                    is_in_core[neighbor] = False
                    q.append(neighbor)
        return is_in_core

    @staticmethod
    def _count_components(num_nodes, adj_list, is_in_core):
        """
        Counts the number of connected components in the resulting k-core using BFS.
        """
        count = 0
        visited = [False] * num_nodes
        for i in range(num_nodes):
            if is_in_core[i] and not visited[i]:
                count += 1
                q = deque([i])
                visited[i] = True
                while q:
                    u = q.popleft()
                    for neighbor in adj_list[u]:
                        if is_in_core[neighbor] and not visited[neighbor]:
                            visited[neighbor] = True
                            q.append(neighbor)
        return count

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            num_neighbors = len(neighbors)

            if num_neighbors < k:
                continue
            
            # --- Step 1: Build the neighbor-induced subgraph G[N(v)] ---
            # This robust implementation explicitly tracks added edges to prevent
            # any possibility of double-counting, which was the root cause of
            # the persistent bug with the test cases.
            node_map = {node_id: i for i, node_id in enumerate(neighbors)}
            
            subgraph_adj = [[] for _ in range(num_neighbors)]
            subgraph_degrees = [0] * num_neighbors
            added_edges = set()

            # Iterate through each node in the subgraph-to-be
            for u_original, u_local in node_map.items():
                # Iterate through its neighbors in the original graph G
                for w_original in G.adj_list[u_original]:
                    # If the neighbor is also part of the subgraph...
                    if w_original in node_map:
                        w_local = node_map[w_original]
                        
                        # To avoid double-counting, represent the edge in a canonical form
                        # (sorted tuple of local IDs) and check if it's already been added.
                        edge = tuple(sorted((u_local, w_local)))
                        
                        if u_local != w_local and edge not in added_edges:
                            subgraph_adj[u_local].append(w_local)
                            subgraph_adj[w_local].append(u_local)
                            subgraph_degrees[u_local] += 1
                            subgraph_degrees[w_local] += 1
                            added_edges.add(edge)

            # --- Step 2: Find k-core using peeling algorithm ---
            core_nodes = kCoreBaseStructuralDiversity._peel_subgraph(
                num_neighbors, subgraph_adj, list(subgraph_degrees), k
            )
            
            # --- Step 3: Count connected components within k-core ---
            sd[v] = kCoreBaseStructuralDiversity._count_components(
                num_neighbors, subgraph_adj, core_nodes
            )

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
