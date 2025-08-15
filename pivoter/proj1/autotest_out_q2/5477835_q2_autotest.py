#!/usr/bin/env python3
# Auto-generated for 5477835

STUDENT_ID = "5477835"
STUDENT_NAME = "Robert Ng"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n


        for u, edges in enumerate(G.adj_list):
            # Grab neighbour graph
            in_neighbour_subgraph = [False] * n
            for v in edges:
                in_neighbour_subgraph[v] = True

            # Calculate community degrees
            community_degrees = [0] * n
            for v in edges:
                for w in G.adj_list[v]:
                    if w != u and in_neighbour_subgraph[w]:
                        community_degrees[w] += 1

            # Peel to get k-core
            removed = [False] * n
            to_peel = deque()
            for v in edges:
                if community_degrees[v] < k:
                    removed[v] = True
                    to_peel.append(v)

            while to_peel:
                v = to_peel.popleft()
                for w in G.adj_list[v]:
                    if w != u and in_neighbour_subgraph[w]:
                        community_degrees[w] -= 1
                        if community_degrees[w] < k and not removed[w]:
                            removed[w] = True
                            to_peel.append(w)


            # BFS from everywhere to get number of components
            sd[u] = 0

            visited = [False] * n
            component_queue = deque()
            for v in edges:
                if visited[v] or removed[v] or not in_neighbour_subgraph[v]:
                    continue
                component_queue.append(v)
                visited[v] = True
                sd[u] += 1

                while component_queue:
                    w = component_queue.popleft()
                    for x in G.adj_list[w]:
                        if in_neighbour_subgraph[x] and not removed[x] and not visited[x]:
                            visited[x] = True
                            component_queue.append(x)


        return sd


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
