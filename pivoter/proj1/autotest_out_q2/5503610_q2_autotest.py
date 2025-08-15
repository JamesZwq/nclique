#!/usr/bin/env python3
# Auto-generated for 5503610

STUDENT_ID = "5503610"
STUDENT_NAME = "Yang Bai"

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
        sd = [0] * n

        for v in range(n):
            neighbour_graph, degree = kCoreBaseStructuralDiversity._construct_neighbour_graph(G, v)
            # empty neighbour graph
            if len(neighbour_graph) == 0:
                continue
            deleted_nodes = kCoreBaseStructuralDiversity._k_core_peeling(neighbour_graph, degree, k)
            count = kCoreBaseStructuralDiversity._count_connected_components(neighbour_graph, deleted_nodes)
            sd[v] = count
        return sd

    @staticmethod
    def _construct_neighbour_graph(G, v):
        v_neighbours = G.adj_list[v]
        vertex_num = len(v_neighbours)
        adj_list = [[] for _ in range(vertex_num)]
        degree = [0] * vertex_num
        # construct the mapping of node to index
        node_to_index = {node: i for i, node in enumerate(v_neighbours)}
        v_neighbour_set = set(v_neighbours)

        for i, u in enumerate(v_neighbours):
            for w in G.adj_list[u]:
                if w in v_neighbour_set:
                    adj_list[i].append(node_to_index[w])
                    degree[i] += 1
        
        return adj_list, degree

    @staticmethod
    def _k_core_peeling(adj_list, degree, k):
        vertex_num = len(adj_list)
        deleted_nodes = [False] * vertex_num

        # push the nodes with degree less than k into the stack
        queue = deque()
        for i, deg in enumerate(degree):
            if deg < k:
                queue.append(i)

        while queue:
            i = queue.popleft()
            if deleted_nodes[i]:
                continue

            deleted_nodes[i] = True
            for j in adj_list[i]:
                if not deleted_nodes[j]:
                    degree[j] -= 1
                    if degree[j] < k:
                        queue.append(j)

        return deleted_nodes

    @staticmethod
    def _count_connected_components(adj_list, deleted_nodes):
        count = 0
        vertex_num = len(adj_list)
        visited = [False] * vertex_num

        # perform dfs to count the connected components
        # i is not the deleted node
        def dfs(i, adj_list, visited):
            for j in adj_list[i]:
                if not deleted_nodes[j] and not visited[j]:
                    visited[j] = True
                    dfs(j, adj_list, visited)

        for i in range(vertex_num):
            if not visited[i] and not deleted_nodes[i]:
                dfs(i, adj_list, visited)
                count += 1
        
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
