#!/usr/bin/env python3
# Auto-generated for 5578253

STUDENT_ID = "5578253"
STUDENT_NAME = "Fanzhuo Duan"

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
        n = G.vertex_num
        adj = G.adj_list      # adj: adjacent_list
        sd = [0] * n

        # traverse each vertex
        for v in range(n):
            nbrs = adj[v]
            # boundary treatment
            # If the total number of neighbors of node v is less than k, then directly return 0.
            # such as graphs with no edges or isolated vertices will be handled correctly.
            if len(nbrs) < k:
                sd[v] = 0
                continue

            set_nbrs = set(nbrs)

            # 1. calculate the degree of each node in N(v) in the subgraph
            deg = {}
            for node in nbrs:
                deg[node] = 0   # initialize the dictionary(degree=0)
                for neighbor in adj[node]:
                    if neighbor in set_nbrs:
                        deg[node] += 1

            # 2. remove all nodes whose degree < k, updating neighbours' degrees
            queue = deque()
            for node in nbrs:
                if deg[node] < k:
                    queue.append(node)   # put all vertices with degree < k into the queue

            removed_node = set()
            while queue:
                node = queue.popleft()
                if node in removed_node:
                    continue
                removed_node.add(node)     # record the deleted nodes
                for neighbor in adj[node]:
                    if neighbor in set_nbrs and neighbor not in removed_node:
                        deg[neighbor] -= 1     # update the degree of the deleted node's neighbors
                        if deg[neighbor] == k - 1:
                            queue.append(neighbor)   # delete the neighbor node if the updated degree < k

            # 3. Remaining vertices form the k-core
            k_core = []
            for node in nbrs:
                if node not in removed_node:
                    k_core.append(node)      # final k-core subgraph node

            # 4. count connected components of k-core(DFS)
            visited = set()
            comp_num = 0
            for node in k_core:
                if node in visited:
                    continue
                comp_num += 1

                stack = [node]
                visited.add(node)
                while stack:
                    a = stack.pop()
                    for neighbor in adj[a]:
                        if neighbor in set_nbrs and neighbor not in removed_node and neighbor not in visited:
                            visited.add(neighbor)
                            stack.append(neighbor)
            sd[v] = comp_num
        return sd

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
