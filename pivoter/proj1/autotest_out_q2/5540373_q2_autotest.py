#!/usr/bin/env python3
# Auto-generated for 5540373

STUDENT_ID = "5540373"
STUDENT_NAME = "Zhuoqi Chen"

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
        def extract_induced_subgraph(G, nodes):
            node_list = sorted(nodes)
            old_to_new = {old: new for new, old in enumerate(node_list)}
            edge_list = []

            for u in node_list:
                for v in G.adj_list[u]:
                    if v in old_to_new and old_to_new[u] < old_to_new[v]:
                        edge_list.append((old_to_new[u], old_to_new[v]))

            return UndirectedUnweightedGraph([(len(node_list), len(edge_list))] + edge_list)

        def k_core_decomposition(G, k):
            n = G.vertex_num
            deg = []
            for i in range(n):
              deg.append(len(G.adj_list[i]))
            visited = [False] * n
            changed = True
            while changed:
                changed = False
                for u in range(n):
                    if not visited[u] and deg[u] < k:
                        visited[u] = True
                        for v in G.adj_list[u]:
                            deg[v] -= 1
                        changed = True
            return visited

        def count_connected_components(G, valid_nodes):
            visited = set()
            count = 0

            def dfs(u):
                stack = [u]
                while stack:
                    node = stack.pop()
                    if node not in visited:
                        visited.add(node)
                        for nbr in G.adj_list[node]:
                            if nbr in valid_nodes and nbr not in visited:
                                stack.append(nbr)

            for u in valid_nodes:
                if u not in visited:
                    dfs(u)
                    count += 1
            return count

        tau = [0] * G.vertex_num

        for v in range(G.vertex_num):
            nbrs = G.adj_list[v]
            if not nbrs:
                continue
            subG = extract_induced_subgraph(G, nbrs)
            core_mask = k_core_decomposition(subG, k)
            valid_nodes = {u for u in range(subG.vertex_num) if not core_mask[u]}
            tau[v] = count_connected_components(subG, valid_nodes)

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
