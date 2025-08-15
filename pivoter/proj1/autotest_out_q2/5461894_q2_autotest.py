#!/usr/bin/env python3
# Auto-generated for 5461894

STUDENT_ID = "5461894"
STUDENT_NAME = "Jingyan Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque,Counter
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

        if not G or k < 0:
            return []

        # TODO
        n = G.vertex_num
        sd = [0] * n
        
        for vertex in range(n):
            adjacent_nodes = G.adj_list[vertex]
            
            if len(adjacent_nodes) == 0:
                sd[vertex] = 0
                continue

            neighbour_nodes = set(adjacent_nodes)
            
            induced_graph = dict()
            for node in neighbour_nodes:
                induced_graph[node] = []
            
            for src in neighbour_nodes:
                for dst in G.adj_list[src]:
                    if dst in neighbour_nodes:
                        if src not in induced_graph:
                            induced_graph[src] = []
                        induced_graph[src].append(dst)

            result = kCoreBaseStructuralDiversity._count_k_core(induced_graph, k)
            sd[vertex] = result


        return sd

    @staticmethod
    def _count_k_core(sub_G, k):
        if not sub_G:
            return 0

        degrees = kCoreBaseStructuralDiversity._initialize_degrees(sub_G)

        to_remove = kCoreBaseStructuralDiversity.prune_low_degree_vertices(sub_G, degrees, k)

        remaining = [n for n in sub_G if n not in to_remove]
        if not remaining:
            return 0

        return kCoreBaseStructuralDiversity._count_connected_components(sub_G, remaining, to_remove)
    
    @staticmethod
    def _initialize_degrees(graph):
        return {node: len(neighbors) for node, neighbors in graph.items()}

    @staticmethod
    def prune_low_degree_vertices(graph, degree, k_threshold):
        marked_for_removal = set()
        candidates = [v for v in graph if degree[v] < k_threshold]

        while candidates:
            new_candidates = []
            for v in candidates:
                if v in marked_for_removal:
                    continue
                marked_for_removal.add(v)
                for neighbor in graph.get(v, []):
                    if neighbor in marked_for_removal:
                        continue
                    degree[neighbor] -= 1
                    if degree[neighbor] < k_threshold:
                        new_candidates.append(neighbor)
            candidates = new_candidates

        return marked_for_removal

    @staticmethod
    def _count_connected_components(graph, valid_nodes, removed):
        visited = set()
        count = 0

        def bfs(start):
            queue = deque([start])
            visited.add(start)
            while queue:
                node = queue.popleft()
                for neighbor in graph.get(node, []):
                    if neighbor not in visited and neighbor not in removed:
                        visited.add(neighbor)
                        queue.append(neighbor)

        for node in valid_nodes:
            if node not in visited:
                bfs(node)
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
