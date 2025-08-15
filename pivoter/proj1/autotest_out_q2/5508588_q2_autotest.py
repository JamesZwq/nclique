#!/usr/bin/env python3
# Auto-generated for 5508588

STUDENT_ID = "5508588"
STUDENT_NAME = "Ziyu Xiong"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
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
            neighbours = G.adj_list[v]
            if not neighbours:
                sd[v] = 0
                continue

            # Construct subgraphs with neighbors as nodes
            neighbours = set(neighbours)
            sub_graph = defaultdict(list)
            for u in neighbours:
                for i in G.adj_list[u]:
                    if i in neighbours:
                        sub_graph[u].append(i)

            # Calculate the number of connected components of k-cores in this subgraph
            sd[v] = kCoreBaseStructuralDiversity._count_k_core(sub_graph, k)
        return sd

    @staticmethod
    def _count_k_core(sub_graph, k):
        if not sub_graph:
            return 0

        # Delete nodes that do not satisfy the k-core condition
        deleted = kCoreBaseStructuralDiversity._peel_k_core(sub_graph, k)

        remains = [i for i in sub_graph if i not in deleted]
        if not remains:
            return 0

        # Count the number of connected components in the remaining nodes
        return kCoreBaseStructuralDiversity._count_components(sub_graph, deleted)

    @staticmethod
    def _peel_k_core(sub_graph, k):
        # Initialize the degree of each node
        degrees = {u: len(neighbours) for u, neighbours in sub_graph.items()}
        deleted = set()     # Deleted nodes
        vertices = deque()  # Candidate deletion queue

        # Initially join nodes with degree less than k
        for u in sub_graph:
            if degrees[u] < k:
                vertices.append(u)
                deleted.add(u)

        # Loop to remove nodes of degree < k and update degrees of neighboring nodes
        while vertices:
            u = vertices.popleft()
            for i in sub_graph[u]:
                if i not in deleted:
                    degrees[i] -= 1
                    if degrees[i] < k:
                        deleted.add(i)
                        vertices.append(i)

        return deleted

    @staticmethod
    def _count_components(sub_graph, deleted):
        visited = set()
        k_core_count = 0

        # Do BFS on each undeleted node, counting the number of connected blocks
        for node in sub_graph:
            if node not in deleted and node not in visited:
                k_core_count += 1
                queue = deque()
                queue.append(node)
                visited.add(node)
                while queue:
                    i = queue.popleft()
                    for j in sub_graph[i]:
                        if j not in deleted and j not in visited:
                            visited.add(j)
                            queue.append(j)
        return k_core_count   # Returns the number of connected components


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
