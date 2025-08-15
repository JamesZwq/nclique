#!/usr/bin/env python3
# Auto-generated for 5491560

STUDENT_ID = "5491560"
STUDENT_NAME = "Aswin Selvakumar"

# ======= 学生代码 =======
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        def k_core_subgraph(vertices):

            # Build the subgraph induced by given vertices
            index_map = {v: i for i, v in enumerate(vertices)}
            rev_map = {i: v for v, i in index_map.items()}
            n = len(vertices)
            adj = [[] for _ in range(n)]
            for v in vertices:
                for u in G.adj_list[v]:
                    if u in index_map:
                        adj[index_map[v]].append(index_map[u])

            # Peeling process to compute k-core
            deg = [len(adj[i]) for i in range(n)]
            removed = [False] * n
            changed = True
            while changed:
                changed = False
                for i in range(n):
                    if not removed[i] and deg[i] < k:
                        removed[i] = True
                        changed = True
                        for nei in adj[i]:
                            deg[nei] -= 1

            # Collect remaining vertices and build adjacency
            remaining = [i for i in range(n) if not removed[i]]
            visited = [False] * n
            components = []

            def dfs(u, comp):
                visited[u] = True
                comp.append(rev_map[u])
                for v in adj[u]:
                    if not removed[v] and not visited[v]:
                        dfs(v, comp)

            for i in remaining:
                if not visited[i]:
                    comp = []
                    dfs(i, comp)
                    components.append(frozenset(comp))
            return components

        tau = [0] * G.vertex_num
        for v in range(G.vertex_num):
            neighbors = G.adj_list[v]
            if not neighbors:
                tau[v] = 0
                continue
            kcores = k_core_subgraph(neighbors)
            tau[v] = len(kcores)
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
