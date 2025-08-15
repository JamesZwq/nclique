#!/usr/bin/env python3
# Auto-generated for 5596805

STUDENT_ID = "5596805"
STUDENT_NAME = "Linfei Du"

# ======= 学生代码 =======
from collections import deque

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
        n = G.vertex_num
        sd = [0] * n

        # Iterate through each vertex v in the graph
        for v_orig in range(n):
            # For each vertex v_orig, compute its k-core structural diversity
            sd[v_orig] = kCoreBaseStructuralDiversity._compute_diversity_for_vertex(G, v_orig, k)

        return sd
    
    @staticmethod
    def _count_components(num_nodes, adj_list, node_set):
        """
        Count the number of connected components in a given node set using BFS.
        """
        if not node_set:
            return 0

        count = 0
        visited = [False] * num_nodes
        
        for i in node_set:
            if not visited[i]:
                count += 1
                q = deque([i])
                visited[i] = True
                while q:
                    u = q.popleft()
                    for v_neighbor in adj_list[u]:
                        if v_neighbor in node_set and not visited[v_neighbor]:
                            visited[v_neighbor] = True
                            q.append(v_neighbor)
        return count

    @staticmethod
    def _compute_diversity_for_vertex(G, v_orig, k):
        """
        Compute k-core structural diversity for a single vertex.
        """
        # Deduplicate neighbor list to ensure processing unique neighbors
        neighbors = sorted(list(set(G.adj_list[v_orig])))
        num_neighbors = len(neighbors)

        if num_neighbors < k:
            return 0

        # Construct neighbor-induced subgraph
        orig_to_sub = {node_id: i for i, node_id in enumerate(neighbors)}
        
        sub_adj_list = [[] for _ in range(num_neighbors)]
        sub_degrees = [0] * num_neighbors
        
        # Use set for fast lookup to improve efficiency
        neighbors_set = set(neighbors)

        for u_sub, u_orig in enumerate(neighbors):
            for w_orig in G.adj_list[u_orig]:
                if w_orig in neighbors_set:
                    if u_orig < w_orig: # Use original ID to avoid duplicate edge processing
                        w_sub = orig_to_sub[w_orig]
                        sub_adj_list[u_sub].append(w_sub)
                        sub_adj_list[w_sub].append(u_sub)
                        sub_degrees[u_sub] += 1
                        sub_degrees[w_sub] += 1
        
        # Compute k-core of the subgraph
        if k == 0:
            return kCoreBaseStructuralDiversity._count_components(num_neighbors, sub_adj_list, set(range(num_neighbors)))

        current_degrees = list(sub_degrees)
        removed = [False] * num_neighbors
        queue = deque()

        for i in range(num_neighbors):
            if current_degrees[i] < k:
                queue.append(i)
                removed[i] = True

        while queue:
            u = queue.popleft()
            for v_neighbor in sub_adj_list[u]:
                if not removed[v_neighbor]:
                    current_degrees[v_neighbor] -= 1
                    if current_degrees[v_neighbor] < k:
                        queue.append(v_neighbor)
                        removed[v_neighbor] = True
        
        k_core_nodes = {i for i, is_removed in enumerate(removed) if not is_removed}

        if not k_core_nodes:
            return 0
            
        # Count connected components in the k-core
        return kCoreBaseStructuralDiversity._count_components(num_neighbors, sub_adj_list, k_core_nodes)

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
