#!/usr/bin/env python3
# Auto-generated for 5459644

STUDENT_ID = "5459644"
STUDENT_NAME = "Yize Shen"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque,defaultdict
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
        tau = [0] * n

        # Precompute degrees for all vertices
        degrees = [len(adj) for adj in G.adj_list]

        for v in range(n):
            # Get neighbors of v
            neighbors = G.adj_list[v]

            if not neighbors:
                tau[v] = 0
                continue

            # Create subgraph induced by neighbors of v
            neighbor_set = set(neighbors)
            subgraph = {}
            degrees_subgraph = {}

            # Build the subgraph and degrees
            for u in neighbors:
                subgraph[u] = [x for x in G.adj_list[u] if x in neighbor_set]
                degrees_subgraph[u] = len(subgraph[u])

            # Compute k-core of the neighbor subgraph
            # Using the standard k-core decomposition algorithm
            queue = deque()
            deleted = set()

            # Initialize queue with vertices having degree < k
            for u in degrees_subgraph:
                if degrees_subgraph[u] < k:
                    queue.append(u)

            # Process the queue
            while queue:
                u = queue.popleft()
                if u in deleted:
                    continue
                deleted.add(u)

                for w in subgraph[u]:
                    if w not in deleted:
                        degrees_subgraph[w] -= 1
                        if degrees_subgraph[w] < k:
                            queue.append(w)

            # The remaining vertices form the k-core
            k_core = [u for u in degrees_subgraph if u not in deleted]

            if not k_core:
                tau[v] = 0
                continue

            # find connected components in the k-core
            visited = set()
            component_count = 0

            for u in k_core:
                if u not in visited:
                    # BFS to find the connected component
                    queue = deque([u])
                    visited.add(u)
                    component_count += 1

                    while queue:
                        current = queue.popleft()
                        for w in subgraph[current]:
                            if w in k_core and w not in visited:
                                visited.add(w)
                                queue.append(w)

            tau[v] = component_count

        return tau

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
