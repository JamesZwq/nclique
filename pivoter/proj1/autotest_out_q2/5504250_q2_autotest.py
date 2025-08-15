#!/usr/bin/env python3
# Auto-generated for 5504250

STUDENT_ID = "5504250"
STUDENT_NAME = "Lance Wang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
# ... existing code ...
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
        n = G.vertex_num
        sd = [0] * n
        for i in range(n):
            sd[i] = kCoreBaseStructuralDiversity._compute_tau_k(G, k, i)
        return sd


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def _compute_tau_k(G, k, v_id):
        neighbors = G.adj_list[v_id]
        if not neighbors:
            return 0

        neighbor_set = set(neighbors)
        node_map = {node: i for i, node in enumerate(neighbors)}
        num_subgraph_nodes = len(neighbors)

        subgraph_adj = [[] for _ in range(num_subgraph_nodes)]
        subgraph_degrees = [0] * num_subgraph_nodes

        for i, u_orig in enumerate(neighbors):
            for v_orig in G.adj_list[u_orig]:
                if v_orig in neighbor_set:
                    # To avoid adding edges twice and self-loops
                    if u_orig < v_orig:
                        j = node_map[v_orig]
                        subgraph_adj[i].append(j)
                        subgraph_adj[j].append(i)
                        subgraph_degrees[i] += 1
                        subgraph_degrees[j] += 1

        # K-core decomposition
        queue = deque([i for i, deg in enumerate(subgraph_degrees) if deg < k])
        removed = [False] * num_subgraph_nodes
        for i in queue:
            removed[i] = True
        
        while queue:
            u = queue.popleft()
            for v_neighbor in subgraph_adj[u]:
                if not removed[v_neighbor]:
                    subgraph_degrees[v_neighbor] -= 1
                    if subgraph_degrees[v_neighbor] < k:
                        queue.append(v_neighbor)
                        removed[v_neighbor] = True
        
        # Count connected components in the remaining graph
        num_components = 0
        visited = [False] * num_subgraph_nodes
        for i in range(num_subgraph_nodes):
            if not removed[i] and not visited[i]:
                num_components += 1
                q_comp = deque([i])
                visited[i] = True
                while q_comp:
                    u_comp = q_comp.popleft()
                    for v_comp in subgraph_adj[u_comp]:
                        if not removed[v_comp] and not visited[v_comp]:
                            visited[v_comp] = True
                            q_comp.append(v_comp)
        
        return num_components

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
