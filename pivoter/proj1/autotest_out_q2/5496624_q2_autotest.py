#!/usr/bin/env python3
# Auto-generated for 5496624

STUDENT_ID = "5496624"
STUDENT_NAME = "Xinbo Li"

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
        List[int]  # tau_k(v) for all v
        """

        def compute_k_core(adj, k):
            # 构造每个节点的度
            deg = {}
            for v in adj:
                deg[v] = len(adj[v])

            # 小于k的node入队
            q = deque()
            for v in deg:
                if deg[v] < k:
                    q.append(v)

            while len(q) > 0:
                v = q.popleft()
                for u in adj[v]:
                    if deg[u] >= k:
                        deg[u] -= 1
                        if deg[u] == k - 1:
                            q.append(u)
                deg[v] = -1

            # k-core，剩余度 >= k 的节点
            result = set()
            for v in deg:
                if deg[v] >= k:
                    result.add(v)
            return result

        def get_c_c(nodes, adj):
            visited = set()
            components = []

            for node in nodes:
                if node not in visited:
                    stack = [node]
                    comp = []

                    while len(stack) > 0:
                        curr = stack.pop()
                        if curr not in visited:
                            visited.add(curr)
                            comp.append(curr)

                            # 把 curr 的邻居里面也在 nodes里的选出
                            neighbors_nodes = []
                            for n in adj[curr]:
                                if n in nodes:
                                    neighbors_nodes.append(n)

                            for nbr in neighbors_nodes:
                                stack.append(nbr)

                    components.append(comp)

            return components

        n = G.vertex_num
        t = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]


            neighbor_graph = {}
            for u in neighbors:
                neighbor_graph[u] = set()
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbor_graph:
                        neighbor_graph[u].add(w)

            # 计算 k-core 中的节点
            kcore_nodes = compute_k_core(neighbor_graph, k)
            components = get_c_c(kcore_nodes, neighbor_graph)
            t[v] = len(components)

        return t

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
