#!/usr/bin/env python3
# Auto-generated for 5550309

STUDENT_ID = "5550309"
STUDENT_NAME = "Bozhi Meng"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def calc_neighbour_k_core(G, center_node, k):
        # get the neighbours
        neighbours = set(G.adj_list[center_node])
        degrees = [0] * G.vertex_num
        # O(|V|)
        for node in neighbours:
            degrees[node] = len([n for n in G.adj_list[node] if n in neighbours])
        # Select the nodes which degree less than k
        # O(|V|)
        unuse_nodes = [node for node in neighbours if degrees[node] < k]
        # O(|V| + |E|)
        while unuse_nodes:
            node = unuse_nodes.pop()
            degrees[node] = 0
            if node in neighbours:
                # Remove from neighbours
                neighbours.remove(node)
            for neighbour in G.adj_list[node]:
                if degrees[neighbour] > 0:
                    degrees[neighbour] -= 1
                    if degrees[neighbour] == k - 1:
                        degrees[neighbour] = 0
                        unuse_nodes.append(neighbour)
        visit = [False] * G.vertex_num
        k_core_count = 0
        # Count k-core in the rest neighbours
        # O(V + E)
        for node in neighbours:
            if visit[node]:
                continue
            k_core_count += 1
            q = deque()
            q.append(node)
            while q:
                u = q.pop()
                if visit[u]:
                    continue
                visit[u] = True
                for v in G.adj_list[u]:
                    if v in neighbours:
                        q.append(v)
        return k_core_count

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        # Time complexity: n loops, calc k-core O(V+E), so the time complexity is O(V(V+E))
        for i in range(n):
            sd[i] = kCoreBaseStructuralDiversity.calc_neighbour_k_core(G, i, k)
        return sd


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
