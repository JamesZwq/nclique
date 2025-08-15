#!/usr/bin/env python3
# Auto-generated for 5469191

STUDENT_ID = "5469191"
STUDENT_NAME = "Jingyue Hu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute the k-core–based structural diversity τ_k(v) for each vertex v in G.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
            The core threshold.

        Returns
        -------
        List[int]
            A list of length G.vertex_num where the v-th entry is τ_k(v), the number
            of k-cores in the neighbor-induced subgraph of vertex v.
        """
        n = G.vertex_num
        # This will hold τ_k(v) for each v
        diversity = [0] * n

        # Process each vertex independently
        for v in range(n):
            # Extract the neighbor set N(v)
            neighbor_set = set(G.adj_list[v])
            # If v has no neighbors, its diversity is zero
            if not neighbor_set:
                diversity[v] = 0
                continue

            # Build degree map for the induced subgraph on N(v)
            # deg[u] = number of neighbors of u that are also in neighbor_set
            deg = {u: 0 for u in neighbor_set}
            for u in neighbor_set:
                for w in G.adj_list[u]:
                    if w in neighbor_set:
                        deg[u] += 1

            # Peel away vertices whose degree < k (k-core decomposition)
            queue = deque(u for u in neighbor_set if deg[u] < k)
            while queue:
                u = queue.popleft()
                # Skip if already removed
                if u not in neighbor_set:
                    continue
                neighbor_set.remove(u)
                # Decrease the degree of its neighbors in the subgraph
                for w in G.adj_list[u]:
                    if w in neighbor_set:
                        deg[w] -= 1
                        # If w now falls below the threshold, schedule it for removal
                        if deg[w] < k:
                            queue.append(w)

            # Count connected components in the remaining subgraph
            visited = set()
            components = 0
            for u in neighbor_set:
                if u not in visited:
                    components += 1
                    bfs_queue = deque([u])
                    visited.add(u)
                    # Standard BFS to mark this component
                    while bfs_queue:
                        x = bfs_queue.popleft()
                        for w in G.adj_list[x]:
                            if w in neighbor_set and w not in visited:
                                visited.add(w)
                                bfs_queue.append(w)

            # Record the number of k-core components
            diversity[v] = components

        return diversity

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
