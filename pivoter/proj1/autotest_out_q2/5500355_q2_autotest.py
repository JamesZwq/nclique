#!/usr/bin/env python3
# Auto-generated for 5500355

STUDENT_ID = "5500355"
STUDENT_NAME = "Yuan Wei"

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
            neighbors = G.adj_list[v]

            if not neighbors:
                sd[v] = 0
                continue

            neighbor_set = set(neighbors)
            subgraph_edges = []

            for u in neighbors:
                for neighbor in G.adj_list[u]:
                    if neighbor in neighbor_set and neighbor > u:
                        subgraph_edges.append((u, neighbor))

            if not subgraph_edges:
                sd[v] = 0
                continue

            vertex_map = {u: i for i, u in enumerate(neighbor_set)}
            mapped_edges = [(vertex_map[u], vertex_map[v]) for u, v in subgraph_edges]

            subgraph = UndirectedUnweightedGraph([(len(neighbor_set), len(mapped_edges))] + mapped_edges)
            k_cores = kCoreBaseStructuralDiversity._get_k_cores(subgraph, k)
            sd[v] = len(k_cores)

        return sd

    @staticmethod
    def _get_k_cores(G, k):
        if k <= 0:
            return kCoreBaseStructuralDiversity._get_connected_components(G)

        adj = [set(neighbors) for neighbors in G.adj_list]
        degrees = [len(neighbors) for neighbors in adj]

        queue = [v for v in range(G.vertex_num) if degrees[v] < k and degrees[v] > 0]

        while queue:
            v = queue.pop()
            for u in adj[v]:
                adj[u].remove(v)
                degrees[u] -= 1
                if degrees[u] == k - 1:
                    queue.append(u)
            adj[v].clear()
            degrees[v] = 0

        visited = [False] * G.vertex_num
        components = []

        for v in range(G.vertex_num):
            if degrees[v] >= k and not visited[v]:
                component = set()
                stack = [v]
                visited[v] = True

                while stack:
                    node = stack.pop()
                    component.add(node)
                    for neighbor in adj[node]:
                        if not visited[neighbor]:
                            visited[neighbor] = True
                            stack.append(neighbor)

                if component:
                    components.append(component)

        return components

    @staticmethod
    def _get_connected_components(G):
        visited = [False] * G.vertex_num
        components = []

        for v in range(G.vertex_num):
            if not visited[v] and G.adj_list[v]:
                component = set()
                stack = [v]
                visited[v] = True

                while stack:
                    node = stack.pop()
                    component.add(node)
                    for neighbor in G.adj_list[node]:
                        if not visited[neighbor]:
                            visited[neighbor] = True
                            stack.append(neighbor)

                if component:
                    components.append(component)

        return components
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
