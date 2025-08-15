#!/usr/bin/env python3
# Auto-generated for 5529330

STUDENT_ID = "5529330"
STUDENT_NAME = "Haomin Li"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n

        for i in range(n):
          # Get the neighbors of vertex i
            neighbor_list = G.adj_list[i]

            if len(neighbor_list) == 0:
                continue

            # Convert to set for faster lookup
            neighbor_set = set(neighbor_list)
            sub_G = {}
            # Build the adjacency list for the induced subgraph
            for u in neighbor_set:
                connections = []
                for v in G.adj_list[u]:
                    if v in neighbor_set:
                        connections.append(v)
                sub_G[u] = connections

            # Compute the number of k_core components in the induced subgraph
            sd[i] = kCoreBaseStructuralDiversity._count_k_core(sub_G, k)

        return sd


    @staticmethod
    def _count_k_core(sub_G, k):
        if not sub_G:
            return 0

        removed = set()
        deg_map = {}

        # Degree initialization + initial filtering
        pending = deque()
        for u, neighbors in sub_G.items():
            deg = len(neighbors)
            deg_map[u] = deg
            if deg < k:
                removed.add(u)
                pending.append(u)

        # Perform degree-based pruning
        def prune(pending, removed, deg_map):
            while pending:
                node = pending.popleft()
                for nbr in sub_G[node]:
                    if nbr in removed:
                        continue
                    deg_map[nbr] -= 1
                    if deg_map[nbr] < k:
                        removed.add(nbr)
                        pending.append(nbr)

        prune(pending, removed, deg_map)


        # Collect remaining valid nodes
        active = []
        for node in sub_G:
            if node not in removed:
                active.append(node)

        if not active:
            return 0

        # Count connected components among remaining nodes
        def bfs(start, visited, removed):
            queue = deque([start])
            visited.add(start)
            while queue:
                curr = queue.popleft()
                for nbr in sub_G[curr]:
                    if nbr not in visited and nbr not in removed:
                        visited.add(nbr)
                        queue.append(nbr)
        visited = set()
        count = 0
        for node in active:
            if node not in visited:
                count += 1
                bfs(node, visited, removed)
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
