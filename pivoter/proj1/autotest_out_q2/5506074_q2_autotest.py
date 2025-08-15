#!/usr/bin/env python3
# Auto-generated for 5506074

STUDENT_ID = "5506074"
STUDENT_NAME = "Xiaoyu Liang"

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
            # Get neighbors of vertex v
            neighbors = set(G.adj_list[v])

            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # Build neighbor-induced subgraph
            # Create mapping from original vertex id to new id in subgraph
            neighbor_list = list(neighbors)
            vertex_map = {orig_id: new_id for new_id, orig_id in enumerate(neighbor_list)}

            # Build adjacency list for neighbor-induced subgraph
            subgraph_adj = [[] for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list):
                for neighbor in G.adj_list[u]:
                    if neighbor in vertex_map:
                        subgraph_adj[i].append(vertex_map[neighbor])

            # Find k-cores in the neighbor-induced subgraph
            sd[v] = kCoreBaseStructuralDiversity._count_k_cores(subgraph_adj, k)

        return sd

    @staticmethod
    def _count_k_cores(adj_list, k):
        n = len(adj_list)
        if n == 0:
            return 0

        # Initialize degrees
        degrees = [len(adj_list[i]) for i in range(n)]

        # Remove vertices with degree < k iteratively
        removed = [False] * n
        changed = True

        while changed:
            changed = False
            for v in range(n):
                if not removed[v] and degrees[v] < k:
                    removed[v] = True
                    changed = True
                    # Update degrees of neighbors
                    for neighbor in adj_list[v]:
                        if not removed[neighbor]:
                            degrees[neighbor] -= 1

        # Find connected components in the remaining graph (k-core)
        visited = [False] * n
        num_components = 0

        for v in range(n):
            if not removed[v] and not visited[v]:
                # Start BFS/DFS from this vertex
                num_components += 1
                queue = deque([v])
                visited[v] = True

                while queue:
                    current = queue.popleft()
                    for neighbor in adj_list[current]:
                        if not removed[neighbor] and not visited[neighbor]:
                            visited[neighbor] = True
                            queue.append(neighbor)

        return num_components
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
