#!/usr/bin/env python3
# Auto-generated for 5501721

STUDENT_ID = "5501721"
STUDENT_NAME = "Aijia Zhang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity:
    @staticmethod
    def process(G, k):
        #Calculate structural diversity for each node based on k-core of its 1-hop neighborhood.
        n = G.vertex_num
        adj = G.adj_list
        result = [0] * n

        def build_induced_subgraph(vertices):
            # Create subgraph with edges only among given vertices
            sub_adj = defaultdict(set)
            for u in vertices:
                for v in adj[u]:
                    if v in vertices:
                        sub_adj[u].add(v)
            return sub_adj

        def k_core_decomposition(sub_adj, k):
            # Iteratively remove nodes with degree < k
            deg = {u: len(neigh) for u, neigh in sub_adj.items()}
            q = deque([u for u in sub_adj if deg[u] < k])
            while q:
                u = q.popleft()
                for v in sub_adj[u]:
                    if deg[v] >= k:
                        deg[v] -= 1
                        if deg[v] == k - 1:
                            q.append(v)
                deg[u] = -1
            return set(u for u, d in deg.items() if d >= k)

        def count_components(sub_adj, valid_nodes):
            # Count connected components in the subgraph
            visited = set()
            components = 0

            def dfs(u):
                stack = [u]
                while stack:
                    node = stack.pop()
                    if node not in visited:
                        visited.add(node)
                        stack.extend([v for v in sub_adj[node] if v in valid_nodes and v not in visited])

            for node in valid_nodes:
                if node not in visited:
                    components += 1
                    dfs(node)
            return components

        for v in range(n):
            N_v = set(adj[v])
            if not N_v:
                result[v] = 0
                continue
            sub_adj = build_induced_subgraph(N_v)
            k_core_nodes = k_core_decomposition(sub_adj, k)
            if not k_core_nodes:
                result[v] = 0
            else:
                result[v] = count_components(sub_adj, k_core_nodes)

        return result

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
