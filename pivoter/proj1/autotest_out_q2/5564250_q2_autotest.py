#!/usr/bin/env python3
# Auto-generated for 5564250

STUDENT_ID = "5564250"
STUDENT_NAME = "Huiling Fan"

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

        for vertex in range(n):
            neighbors = set(G.adj_list[vertex])

            if len(neighbors) == 0:
                sd[vertex] = 0
                continue

            neighbors_graph = {}

            for node in neighbors:
                neighbors_graph[node] = []
                for u in G.adj_list[node]:
                    if u in neighbors:
                        neighbors_graph[node].append(u)

            k_core = kCoreBaseStructuralDiversity._find_k_cores(neighbors_graph, k)
            sd[vertex] = k_core

        return sd

    @staticmethod
    def _get_connected_components(graph):
        all_visited = set()
        comp_list = []

        for source in graph:
            if source in all_visited:
                continue

            frontier = deque([source])
            block = {source}
            all_visited.add(source)

            while frontier:
                top = frontier.popleft()
                for adj in graph[top]:
                    if adj not in all_visited:
                        frontier.append(adj)
                        block.add(adj)
                        all_visited.add(adj)

            comp_list.append(block)

        return comp_list


    @staticmethod
    def _count_degv(component, node_subset):
        degrees = {}
        for node in node_subset:
            count = 0
            for neigh in component[node]:
                if neigh in node_subset:
                    count += 1
            degrees[node] = count
        return degrees

    @staticmethod
    def _find_k_cores(graph, k):
        blocks = kCoreBaseStructuralDiversity._get_connected_components(graph)
        count = 0

        for block in blocks:
            
            if len(block) < k:
                continue

           
            subgraph = {}
            for node in block:
                subgraph[node] = [nbr for nbr in graph[node] if nbr in block]

            alive = set(subgraph.keys())

            while True:
                degrees = kCoreBaseStructuralDiversity._count_degv(subgraph, alive)
                to_remove = [node for node in alive if degrees[node] < k]
                if not to_remove:
                    break
                for node in to_remove:
                    alive.remove(node)

            if alive:
                
                filtered_graph = {u: [v for v in subgraph[u] if v in alive] for u in alive}
                subblocks = kCoreBaseStructuralDiversity._get_connected_components(filtered_graph)
                count += len(subblocks)

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
