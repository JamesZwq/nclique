#!/usr/bin/env python3
# Auto-generated for 5461960

STUDENT_ID = "5461960"
STUDENT_NAME = "Yiran Li"

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
        adj = G.adj_list
        tau = [0] * n

        adj_sets = [set(neighbors) for neighbors in adj]
        degrees_list = [len(s) for s in adj_sets]

        adj_directed = [[] for _ in range(n)]
        for u in range(n):
            for v in adj[u]:
                if degrees_list[u] < degrees_list[v] or (degrees_list[u] == degrees_list[v] and u < v):
                    adj_directed[u].append(v)

        triangle_counts = {}
        for u in range(n):
            for i in range(len(adj_directed[u])):
                for j in range(i + 1, len(adj_directed[u])):
                    v = adj_directed[u][i]
                    w = adj_directed[u][j]
                    if w in adj_sets[v]:
                        uv_key = (u, v) if u < v else (v, u)
                        uw_key = (u, w) if u < w else (w, u)
                        vw_key = (v, w) if v < w else (w, v)

                        triangle_counts[uv_key] = triangle_counts.get(uv_key, 0) + 1
                        triangle_counts[uw_key] = triangle_counts.get(uw_key, 0) + 1
                        triangle_counts[vw_key] = triangle_counts.get(vw_key, 0) + 1

        for v in range(n):
            neighbors = adj_sets[v]

            if len(neighbors) < k or k == 0:
                tau[v] = 0
                continue

            degrees_in_subgraph = {u: triangle_counts.get((min(v, u), max(v, u)), 0) for u in neighbors}

            to_remove = deque([u for u, d in degrees_in_subgraph.items() if d < k])
            removed_nodes = set(to_remove)

            while to_remove:
                u_removed = to_remove.popleft()
                for neighbor_of_u in adj[u_removed]:
                    if neighbor_of_u in degrees_in_subgraph and neighbor_of_u not in removed_nodes:
                        degrees_in_subgraph[neighbor_of_u] -= 1
                        if degrees_in_subgraph[neighbor_of_u] == k - 1:
                            removed_nodes.add(neighbor_of_u)
                            to_remove.append(neighbor_of_u)

            core_nodes = neighbors - removed_nodes

            component_count = 0
            visited = set()

            for node in core_nodes:
                if node not in visited:
                    component_count += 1
                    q_bfs = deque([node])
                    visited.add(node)
                    while q_bfs:
                        curr = q_bfs.popleft()
                        for neighbor_of_curr in adj[curr]:
                            if neighbor_of_curr in core_nodes and neighbor_of_curr not in visited:
                                visited.add(neighbor_of_curr)
                                q_bfs.append(neighbor_of_curr)

            tau[v] = component_count

        return tau

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
