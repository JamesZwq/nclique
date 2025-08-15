#!/usr/bin/env python3
# Auto-generated for 5466935

STUDENT_ID = "5466935"
STUDENT_NAME = "Xinyi Chen"

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

        # get neighbor and subgraph
        for v in range(n):
            neighbors = G.adj_list[v]
            if not neighbors:
              sd[v] = 0
              continue

            subgraph_nodes = set(neighbors)

            # get new index for subgraph
            node_index = {}
            i = 0
            for node in subgraph_nodes:
                node_index[node] = i
                i += 1

            # build subgraph adjacency list
            subgraph_adj = []
            for _ in range(len(subgraph_nodes)):
                subgraph_adj.append([])

            for u in subgraph_nodes:
                for v2 in G.adj_list[u]:
                    if v2 in subgraph_nodes:
                        subgraph_adj[node_index[u]].append(node_index[v2])

            # remove nodes degrees < k from the subgraph
            deg = [len(adj) for adj in subgraph_adj]
            removed = [False] * len(subgraph_adj)
            queue = deque()
            for i in range(len(deg)):
                if deg[i] < k:
                    queue.append(i)

            while queue:
                u = queue.popleft()
                if removed[u]:
                    continue
                removed[u] = True
                for v2 in subgraph_adj[u]:
                    if not removed[v2]:
                        deg[v2] -= 1
                        if deg[v2] < k:
                            queue.append(v2)
            if all(removed):
              sd[v] = 0
              continue

            # count connected components (BFS)
            visited = [False] * len(subgraph_adj)
            count = 0
            for i in range(len(subgraph_adj)):
                if not removed[i] and not visited[i]:
                    count += 1
                    q = deque([i])
                    visited[i] = True
                    while q:
                        u = q.popleft()
                        for v2 in subgraph_adj[u]:
                            if not removed[v2] and not visited[v2]:
                                visited[v2] = True
                                q.append(v2)
            sd[v] = count
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
