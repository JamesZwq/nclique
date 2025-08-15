#!/usr/bin/env python3
# Auto-generated for 5516805

STUDENT_ID = "5516805"
STUDENT_NAME = "Haitao Wu"

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
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
                continue

            # build neibor-incloude subgraph using adj_list
            # only includes edges between neibors of v
            subgraph = defaultdict(list)
            neighbor_set = set(neighbors)
            for u in neighbor_set:
                for n in G.adj_list[u]:
                    if n in neighbor_set:
                        subgraph[u].append(n)
            #get the result
            sd[v] = kCoreBaseStructuralDiversity.count_k_core_components(subgraph, k)

        return sd

    @staticmethod
    def count_k_core_components(subgraph, k):
        # compute degree of each node in the subgraph
        degree = {u: len(neis) for u, neis in subgraph.items()}
        # set removal queue with nodes having degree < k
        removed = set()
        queue = deque([u for u in subgraph if degree[u] < k])

        # remove low-degree nodes and update neighbors
        while queue:
            u = queue.popleft()
            if u in removed:
                continue
            removed.add(u)
            for v in subgraph[u]:
                if v not in removed:
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)

        visited = set()
        components = 0

        for u in subgraph:
            # skip the already removed or visited nodes
            if u in removed or u in visited:
                continue
            # use BFS to find all reachable nodes forming a k-core component
            components += 1
            queue = deque([u])
            visited.add(u)
            while queue:
                node = queue.popleft()
                for n in subgraph[node]:
                    if n not in removed and n not in visited:
                        visited.add(n)
                        queue.append(n)

        return components

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
