#!/usr/bin/env python3
# Auto-generated for 5467664

STUDENT_ID = "5467664"
STUDENT_NAME = "Yuxuan Ren"

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
            neighbours = G.adj_list[v]
            if not neighbours:
                continue

            # Build neighbor subgraph
            sub_G = kCoreBaseStructuralDiversity._build_subgraph(G, neighbours)

            # Calculate k-core count
            sd[v] = kCoreBaseStructuralDiversity._compute_k_core(sub_G, k)

        return sd

    @staticmethod
    def _build_subgraph(G, neighbours):
        # Remove duplicate neighbors
        neighbours_set = set(neighbours)

        # Build adjacency list for subgraph
        sub_G = {}
        for u in neighbours_set:
            sub_G[u] = []
            for v in G.adj_list[u]:
                if v in neighbours_set:
                    sub_G[u].append(v)

        return sub_G

    @staticmethod
    def _compute_k_core(sub_G, k):
        if not sub_G:
            return 0

        # Calculate degree for each vertex
        degrees = {}
        for u in sub_G:
            degrees[u] = len(sub_G[u])

        deleted = set()
        queue = deque()

        # Find all vertices with degree less than k
        for u in sub_G:
            if degrees[u] < k:
                queue.append(u)
                deleted.add(u)

        # Remove vertices with insufficient degree and update neighbors' degrees
        while queue:
            u = queue.popleft()
            for v in sub_G[u]:
                if v not in deleted:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        deleted.add(v)
                        queue.append(v)

        # Find remaining vertices
        remaining = []
        for u in sub_G:
            if u not in deleted:
                remaining.append(u)

        if not remaining:
            return 0

        # Count connected components
        visited = set()
        count = 0

        for node in remaining:
            if node not in visited:
                count += 1
                # BFS traverse connected component
                bfs_queue = deque()
                bfs_queue.append(node)
                visited.add(node)

                while bfs_queue:
                    current = bfs_queue.popleft()
                    for neighbor in sub_G[current]:
                        if neighbor not in visited and neighbor not in deleted:
                            visited.add(neighbor)
                            bfs_queue.append(neighbor)

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
