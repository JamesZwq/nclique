#!/usr/bin/env python3
# Auto-generated for 5440601

STUDENT_ID = "5440601"
STUDENT_NAME = "Yaolun Luo"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _find_k_cores_and_count_components(subgraph_adj, degrees, k):
        """
        Finds the k-cores of a subgraph and counts the number of connected components.

        Parameters:
        subgraph_adj: Adjacency list of the subgraph.
        degrees: A list of degrees for each node in the subgraph.
        k: The core number.

        Returns:
        The number of connected components in the k-core of the subgraph.
        """

        # Track removed nodes instead of creating a new graph
        queue = deque()
        removed = [False] * len(subgraph_adj)
        
        # Initialize queue with all nodes having degree less than k
        for i in range(len(degrees)):
            if degrees[i] < k:
                queue.append(i)
                removed[i] = True

        # Iteratively remove nodes with degree less than k
        while queue:
            u = queue.popleft()
            for v in subgraph_adj[u]:
                if not removed[v]:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)
                        removed[v] = True
        
        component_count = 0
        
        for i in range(len(subgraph_adj)):
            if not removed[i]:
                component_count += 1
                # Start BFS to find all nodes in this component
                component_queue = deque([i])
                removed[i] = True
                while component_queue:
                    u = component_queue.popleft()
                    for v in subgraph_adj[u]:
                        # If neighbor v is also in the k-core and not visited
                        if not removed[v]:
                            component_queue.append(v)
                            removed[v] = True
                            
        return component_count

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
        n = G.vertex_num
        # sd stands for structural diversity, this will be our final result list τ
        sd = [0] * n

        # Iterate through each vertex in the graph G
        for i in range(n):
            neighbors = G.adj_list[i]
            num_neighbors = len(neighbors)

            # If a node has fewer than k neighbors then have no k-core can exist.
            if num_neighbors <= k:
                sd[i] = 0
                continue

            # Build the neighbor-induced subgraph
            node_map = {node_id: new_id for new_id, node_id in enumerate(neighbors)}
            subgraph_adj = [[] for _ in range(num_neighbors)]
            subgraph_deg = [0] * num_neighbors
            
            # Iterate through each neighbor u of i
            for u_original in neighbors:
                u_local = node_map[u_original]
                # For each neighbor v of u
                for v_original in G.adj_list[u_original]:
                    # If v is also a neighbor of i, then the edge (u,v) exists in the neighbor-induced subgraph.
                    if v_original in node_map:
                        v_local = node_map[v_original]
                        # To avoid double counting edges and degrees
                        if u_original < v_original:
                            subgraph_adj[u_local].append(v_local)
                            subgraph_adj[v_local].append(u_local)
                            subgraph_deg[u_local] += 1
                            subgraph_deg[v_local] += 1
            
            # Find k-cores and count components
            sd[i] = kCoreBaseStructuralDiversity._find_k_cores_and_count_components(subgraph_adj, list(subgraph_deg), k)
            
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
