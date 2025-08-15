#!/usr/bin/env python3
# Auto-generated for 5519051

STUDENT_ID = "5519051"
STUDENT_NAME = "Chenxi Liu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque


class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):

        n = G.vertex_num
        result = [0] * n

        for v in range(n):
            nbr_subgraph = kCoreBaseStructuralDiversity.get_neighbour_induced_subgraph(G, v)
            k_core = kCoreBaseStructuralDiversity.extract_k_core(nbr_subgraph, k)
            result[v] = kCoreBaseStructuralDiversity.count_connected_components(k_core)

        return result

    
    @staticmethod
    def get_neighbour_induced_subgraph(G, v):
        neighbors = G.adj_list[v]
        neighbor_set = set(neighbors)
        
        # subgraph
        subgraph = {}
        for u in neighbor_set:
            subgraph[u] = set()
        
        # add edges
        for u in neighbor_set:
            for w in G.adj_list[u]:
                if w in neighbor_set:
                    subgraph[u].add(w)
                    subgraph[w].add(u)

        # output -> list
        result = {}
        for u in subgraph:
            result[u] = list(subgraph[u])
            
        return result


    @staticmethod
    def extract_k_core(subgraph, k):

        # degree of nodes
        deg = {}
        for u in subgraph:
            deg[u] = len(subgraph[u])

        # queue for degree < k
        queue = deque()
        for u in subgraph:
            if deg[u] < k:
                queue.append(u)

        # bfs
        while queue:
            u = queue.popleft()
            for v in subgraph[u]:
                if v in deg and deg[v] > 0:
                    deg[v] -= 1
                    if deg[v] == k - 1:
                        queue.append(v)
                    
            deg[u] = 0

        # construct k core subgraph

        core = {}
        for u in subgraph:
            if deg[u] >= k:
                core[u] = []

        # neighbour of nodes
        for u in core:
            for v in subgraph[u]:
                if v in core:
                    core[u].append(v)

        return core



    @staticmethod
    # connected components numbers
    def count_connected_components(graph):
        visited = set()
        count = 0

        for node in graph:
            if node not in visited:
                count += 1
                queue = deque([node])
                visited.add(node)
                while queue:
                    u = queue.popleft()
                    for v in graph[u]:
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)
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
