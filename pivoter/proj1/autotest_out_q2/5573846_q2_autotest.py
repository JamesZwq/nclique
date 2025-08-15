#!/usr/bin/env python3
# Auto-generated for 5573846

STUDENT_ID = "5573846"
STUDENT_NAME = "Wei Ting Lai"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################
from collections import defaultdict

class kCoreBaseStructuralDiversity:
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
        List[int] : τ_k(v) for all v
        """
        n = G.vertex_num
        adj = G.adj_list
        result = [0] * n
# set up the array to store the vertex
        for v in range(n):
# input the vertex v exclude v itself
            neighbors = list(set(adj[v]) - {v})
            if not neighbors:
                continue

# map node IDs for induced subgraph
            node_map = {u: i for i, u in enumerate(neighbors)}
            sub_adj = [[] for _ in range(len(neighbors))]
            for i, u in enumerate(neighbors):
                for w in G.adj_list[u]:
                    if w in node_map:
                        sub_adj[i].append(node_map[w])
# compute connected comppnents in the k-core 
# Run local k-core + count connected components
            result[v] = kCoreBaseStructuralDiversity._local_k_core_count(sub_adj, k)

        return result


# Function remove all nodes which is less than k in the graph and counts the number of it among remaining nodes. 
# adj_list : List[List[int]] Adjacency list of the induced subgraph.
    @staticmethod
    def _local_k_core_count(adj_list, k):
        n = len(adj_list)
        deg = [len(adj_list[i]) for i in range(n)]
        removed = [False] * n
        active = [i for i in range(n) if deg[i] < k]
# remove nodes with degree<k and update neighbors' degrees        
        idx = 0
        while idx < len(active):
            u = active[idx]
            removed[u] = True
            for v in adj_list[u]:
                if not removed[v]:
                    deg[v] -= 1
                    if deg[v] == k - 1:
                        active.append(v)
            idx += 1
# use BFS count connected components
        visited = [False] * n
        count = 0
        for i in range(n):
            if removed[i] or visited[i]:
                continue
            count += 1
            stack = [i]
            visited[i] = True
# use DFS to explore the rest
            while stack:
                u = stack.pop()
                for v in adj_list[u]:
                    if not removed[v] and not visited[v]:
                        visited[v] = True
                        stack.append(v)
        return count

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
