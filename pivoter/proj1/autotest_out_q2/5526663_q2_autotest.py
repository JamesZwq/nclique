#!/usr/bin/env python3
# Auto-generated for 5526663

STUDENT_ID = "5526663"
STUDENT_NAME = "Tianyang Wan"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _count_connected_components(adj, is_removed):
        """
        Counts connected components in a graph using BFS, ignoring removed nodes.
        """
        count = 0
        num_nodes = len(adj)
        visited = [False] * num_nodes
        
        for i in range(num_nodes):
            if not is_removed[i] and not visited[i]:
                count += 1
                q = deque([i])
                visited[i] = True
                while q:
                    u = q.popleft()
                    for v_neighbor in adj[u]:
                        if not is_removed[v_neighbor] and not visited[v_neighbor]:
                            visited[v_neighbor] = True
                            q.append(v_neighbor)
        return count

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
        tau_values = [0] * n

        for v in range(n):
            # Ensure the list of neighbors is unique before any processing.
            # This handles cases where G.adj_list[v] might contain duplicates.
            neighbors = sorted(list(set(G.adj_list[v])))
            
            num_neighbors = len(neighbors)

            # if the number of neighbors is less than k+1,
            # it's impossible to form a k-core.
            if num_neighbors < k + 1:
                tau_values[v] = 0
                continue
            
            # --- Subgraph Construction ---
            node_map = {node_id: i for i, node_id in enumerate(neighbors)}
            neighbor_set = set(neighbors)
            
            subgraph_adj = [[] for _ in range(num_neighbors)]
            subgraph_degrees = [0] * num_neighbors

            for i, u_original in enumerate(neighbors):
                for w_original in G.adj_list[u_original]:
                    if u_original < w_original and w_original in neighbor_set:
                        j = node_map[w_original]
                        subgraph_adj[i].append(j)
                        subgraph_adj[j].append(i)
                        subgraph_degrees[i] += 1
                        subgraph_degrees[j] += 1
            
            # --- Peeling Algorithm ---
            degrees = list(subgraph_degrees)
            is_removed = [False] * num_neighbors
            
            queue = deque()
            for i in range(num_neighbors):
                if degrees[i] < k:
                    queue.append(i)
                    is_removed[i] = True
            
            while queue:
                u_sub = queue.popleft()
                for v_sub in subgraph_adj[u_sub]:
                    if not is_removed[v_sub]:
                        degrees[v_sub] -= 1
                        if degrees[v_sub] < k:
                            is_removed[v_sub] = True
                            queue.append(v_sub)
            
            # --- Counting Connected Components ---
            tau_values[v] = kCoreBaseStructuralDiversity._count_connected_components(subgraph_adj, is_removed)

        return tau_values

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
