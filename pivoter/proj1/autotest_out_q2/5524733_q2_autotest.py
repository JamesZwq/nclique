#!/usr/bin/env python3
# Auto-generated for 5524733

STUDENT_ID = "5524733"
STUDENT_NAME = "Haizheng Miao"

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
        adj = G.adj_list
        result = [0] * n
        
        # For each node v, handle it separately
        for v in range(n):
            neighbors = adj[v]
            if not neighbors:
                result[v] = 0
                continue

            # Step 1. Establish a neighbor induction subgraph
            nbr_set = set(neighbors)
            subgraph_nodes = list(nbr_set)
            node_id_map = {u: i for i, u in enumerate(subgraph_nodes)}  # node -> idx in subgraph
            
            # Construct the adjacency list of the neighbor induction subgraph
            sub_adj = [[] for _ in subgraph_nodes]
            for i, u in enumerate(subgraph_nodes):
                for w in adj[u]:
                    if w in nbr_set:
                        sub_adj[i].append(node_id_map[w])

            # Step 2. K-core stripping (deleting points with a degree less than k until convergence)
            degrees = [len(lst) for lst in sub_adj]
            active = [True] * len(subgraph_nodes)
            queue = deque([i for i, deg in enumerate(degrees) if deg < k])

            while queue:
                u = queue.popleft()
                if not active[u]:
                    continue
                active[u] = False
                for w in sub_adj[u]:
                    if active[w]:
                        degrees[w] -= 1
                        if degrees[w] == k - 1:
                            queue.append(w)

            # Step 3. Perform BFS on the remaining points to find the connected components
            visited = [False] * len(subgraph_nodes)
            count = 0
            for i in range(len(subgraph_nodes)):
                if active[i] and not visited[i]:
                    count += 1
                    dq = deque([i])
                    visited[i] = True
                    while dq:
                        cur = dq.popleft()
                        for nxt in sub_adj[cur]:
                            if active[nxt] and not visited[nxt]:
                                visited[nxt] = True
                                dq.append(nxt)
            result[v] = count
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
