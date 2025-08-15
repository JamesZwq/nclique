#!/usr/bin/env python3
# Auto-generated for 5380003

STUDENT_ID = "5380003"
STUDENT_NAME = "Jia Li"

# ======= 学生代码 =======
from collections import deque, defaultdict

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
        # This function gets the k-core subgraph from the given graph
        def get_k_core(graph_adj):
            degree = {node: len(neigh) for node, neigh in graph_adj.items()}
            queue = deque([node for node in graph_adj if degree[node] < k])

            # Remove nodes that do not satisfy k-core condition
            while queue:
                node = queue.popleft()
                for neigh in graph_adj[node]:
                    if neigh in degree:
                        degree[neigh] -= 1
                        if degree[neigh] == k - 1:
                            queue.append(neigh)
                degree.pop(node)

            kcore_adj = defaultdict(list)
            for node in degree:
                for neigh in graph_adj[node]:
                    if neigh in degree:
                        kcore_adj[node].append(neigh)
            return kcore_adj

        # This function calculates the number of connected components in a graph
        def count_connected_components(adj):
            visited = set()
            count = 0

            for node in adj:
                if node not in visited:
                    count += 1
                    stack = [node]

                    while stack:
                        curr = stack.pop()
                        if curr in visited:
                            continue
                        visited.add(curr)
                        for neigh in adj[curr]:
                            if neigh not in visited:
                                stack.append(neigh)
            return count
        τ = [0] * G.vertex_num

        # Calculate τ_k(v) for each node v in the graph
        for v in range(G.vertex_num):
            neighbors = G.adj_list[v]
            subgraph_adj = defaultdict(list)
            neighbor_set = set(neighbors)
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbor_set:
                        subgraph_adj[u].append(w)

            kcore_adj = get_k_core(subgraph_adj)

            τ[v] = count_connected_components(kcore_adj)

        return τ

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
