#!/usr/bin/env python3
# Auto-generated for 5520557

STUDENT_ID = "5520557"
STUDENT_NAME = "Pengju Chen"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from collections import defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [0] * n
        for v in range(n):
            #build neighbor_induced_subgraph
            neighbors = set(G.adj_list[v])
            subgraph = defaultdict(set)
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbors:
                        subgraph[u].add(w)
                        subgraph[w].add(u)
            #use function k_core to calculate k-core of subgraph
            #subgraph:{u:[neighbors]}
            k_core_subgraph = kCoreBaseStructuralDiversity.get_k_core(subgraph, k)
            #count connect components
            sd[v] = kCoreBaseStructuralDiversity.count_components(k_core_subgraph)
        return sd

    @staticmethod
    def get_k_core(subgraph, k):
        #get degree of each node, degrees:{node:degree}
        degrees = {node: len(neighbors) for node, neighbors in subgraph.items()}
        #add node whose degree is less than k into deque
        queue = deque([node for node in degrees if degrees[node] < k])
        while queue:
            u = queue.popleft()
            for v in subgraph[u]:
                #remove u from v's neighbors list
                subgraph[v].discard(u)
                if len(subgraph[v]) == k - 1:
                    queue.append(v)
            #remove all neighbors of u
            subgraph[u].clear()

        return {u: neighbors for u, neighbors in subgraph.items() if neighbors}

    @staticmethod
    def count_components(graph):
        visited = set()
        count = 0
        #use bfs to traverse
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
