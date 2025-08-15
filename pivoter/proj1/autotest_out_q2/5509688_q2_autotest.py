#!/usr/bin/env python3
# Auto-generated for 5509688

STUDENT_ID = "5509688"
STUDENT_NAME = "Ganming Luo"

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
        # TODOv
        n = G.vertex_num
        sd = [0] * n

        # For each vertex v, compute the number of k-cores in its neighbor-induced subgraph
        for v in range(n):
            neighbours = G.adj_list[v]
            if not neighbours:
                sd[v] = 0
                continue

            # Construct neighbor-induced subgraph
            neighbours_set = set(neighbours)
            sub_G = {u: [] for u in neighbours_set}
            for u in neighbours_set:
                # Only keep neighbors that are in neighbours_set
                sub_G[u] = [w for w in G.adj_list[u] if w in neighbours_set]

            # Calculate the number of k-cores in this subgraph
            sd[v] = kCoreBaseStructuralDiversity._count_k_core(sub_G, k)

        return sd

    @staticmethod
    def _count_k_core(sub_G, k):
        # If subgraph is empty, there are no k-cores
        if not sub_G:
            return 0

        # Initialize degree table and deletion queue
        degrees = {u: len(neigh) for u, neigh in sub_G.items()}
        deleted = set()
        queue = deque()

        # First round: delete nodes with degree < k
        for u, deg in degrees.items():
            if deg < k:
                deleted.add(u)
                queue.append(u)

        # Iteratively delete and update neighbor degrees
        while queue:
            u = queue.popleft()
            for v in sub_G[u]:
                if v not in deleted:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        deleted.add(v)
                        queue.append(v)

        # Remaining nodes are all k-core nodes
        remains = [u for u in sub_G if u not in deleted]
        if not remains:
            return 0

        # Count the number of connected components
        visited = set()
        count = 0
        for u in remains:
            if u not in visited:
                count += 1
                bfs = deque([u])
                visited.add(u)
                while bfs:
                    x = bfs.popleft()
                    for y in sub_G[x]:
                        if y not in deleted and y not in visited:
                            visited.add(y)
                            bfs.append(y)
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
