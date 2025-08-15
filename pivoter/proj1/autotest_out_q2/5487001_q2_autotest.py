#!/usr/bin/env python3
# Auto-generated for 5487001

STUDENT_ID = "5487001"
STUDENT_NAME = "Yekun Wang"

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
        n = G.vertex_num
        result = [0] * n

        for node in range(n):
            neighbors = G.adj_list[node]
            if not neighbors:
                result[node] = 0
                continue

            subgraph = kCoreBaseStructuralDiversity._induced_subgraph(G, neighbors)


            result[node] = kCoreBaseStructuralDiversity._count_k_core_components(subgraph, k)

        return result

    @staticmethod
    def _induced_subgraph(G, node_list):

        sub_adj = dict()
        node_set = set(node_list)
        for u in node_list:
            sub_adj[u] = [v for v in G.adj_list[u] if v in node_set]
        return sub_adj

    @staticmethod
    def _count_k_core_components(sub_adj, k):

        if not sub_adj:
            return 0

        degrees = {u: len(neis) for u, neis in sub_adj.items()}
        to_remove = deque([u for u in sub_adj if degrees[u] < k])
        removed = set(to_remove)

        while to_remove:
            u = to_remove.popleft()
            for v in sub_adj[u]:
                if v not in removed:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        removed.add(v)
                        to_remove.append(v)


        remaining_nodes = [u for u in sub_adj if u not in removed]
        visited = set()
        count = 0

        for u in remaining_nodes:
            if u not in visited:
                count += 1
                queue = deque([u])
                visited.add(u)

                while queue:
                    curr = queue.popleft()
                    for v in sub_adj[curr]:
                        if v not in removed and v not in visited:
                            visited.add(v)
                            queue.append(v)

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
