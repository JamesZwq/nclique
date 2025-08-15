#!/usr/bin/env python3
# Auto-generated for 5408287

STUDENT_ID = "5408287"
STUDENT_NAME = "Geng Wu"

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
        core_numbers = kCoreBaseStructuralDiversity._calculate_core_numbers(G)

        for u in range(n):
            sub_graph = kCoreBaseStructuralDiversity._obtain_subgraph(G, u, core_numbers, k)
            if sub_graph.vertex_num == 0 or sub_graph.vertex_num < k:
                continue
            
            remaining_vertexes = kCoreBaseStructuralDiversity._k_core_computation(sub_graph, k)
            num_k_cores = kCoreBaseStructuralDiversity._obtain_num_k_cores(sub_graph, remaining_vertexes)
            sd[u] = num_k_cores
        return sd
    
    @staticmethod
    def _calculate_core_numbers(G):
        n = G.vertex_num
        degree = [len(G.adj_list[v]) for v in range(n)]

        max_deg = max(degree) if n else 0
        bin_boundaries = [0] * (max_deg + 2)

        for d in degree:
            bin_boundaries[d] += 1
        
        start = 0
        for d in range(max_deg + 1):
            cnt = bin_boundaries[d]
            bin_boundaries[d] = start
            start += cnt

        D = [0] * n
        pos = [0] * n
        next_slot = bin_boundaries[:]

        for v, d in enumerate(degree):
            idx = next_slot[d]
            D[idx] = v
            pos[v] = idx
            next_slot[d] += 1

        for i in range(n):
            v = D[i]
            for u in G.adj_list[v]:
                if degree[u] > degree[v]:
                    du = degree[u]
                    pu = pos[u]
                    pw = bin_boundaries[du]
                    w = D[pw]

                    if u != w:
                        D[pu], D[pw] = D[pw], D[pu]
                        pos[u], pos[w] = pw, pu

                    bin_boundaries[du] += 1
                    degree[u] -= 1

        return degree

    @staticmethod
    def _obtain_subgraph(G, u, core_numbers, k):
        nbr_vertexs = [v for v in G.adj_list[u] if core_numbers[v] >= k]
        new_indexes = {v: i for i, v in enumerate(nbr_vertexs)}
        
        class TargetGraph:
            def __init__(self, nbr_vertexs, new_indexes, G):
                self.vertex_num = len(nbr_vertexs)
                self.adj_list = [list() for _ in range(self.vertex_num)]

                for v1 in nbr_vertexs:
                    self.adj_list[new_indexes[v1]].extend(new_indexes[v2] for v2 in G.adj_list[v1] if v2 in new_indexes)
        
        return TargetGraph(nbr_vertexs, new_indexes, G)

    @staticmethod
    def _k_core_computation(G, k):
        n = G.vertex_num
        degree = [len(G.adj_list[v]) for v in range(n)]

        q = deque(v for v, d in enumerate(degree) if d < k)
        remaining_vertexes = [True] * n

        while len(q) > 0:
            v = q.popleft()
            if not remaining_vertexes[v]:
                continue

            remaining_vertexes[v] = False
            
            for n in G.adj_list[v]:
                if not remaining_vertexes[n]:
                    continue

                degree[n] -= 1
                if degree[n] < k:
                    q.append(n)

        return remaining_vertexes
    
    @staticmethod
    def _obtain_num_k_cores(G, remaining_vertexes):
        
        num_k_cores = 0
        q = deque()
        visited = [False] * G.vertex_num

        for v in range(G.vertex_num):
            if not remaining_vertexes[v]:
                continue

            if visited[v]:
                continue

            num_k_cores += 1
            q.append(v)
            visited[v] = True

            while len(q) > 0:
                u = q.popleft()
                for n in G.adj_list[u]:
                    if not visited[n] and remaining_vertexes[n]:
                        visited[n] = True
                        q.append(n)
        
        return num_k_cores


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
