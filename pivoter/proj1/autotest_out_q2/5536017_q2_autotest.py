#!/usr/bin/env python3
# Auto-generated for 5536017

STUDENT_ID = "5536017"
STUDENT_NAME = "Bonnie Lu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        G: UndirectedUnweightedGraph
        k: int
        Return: List[int]  # τ_k(v) for all v
        """
        n = G.vertex_num
        result = [0] * n

        # Go through each vertex
        for v in range(n):
            nbrs = G.adj_list[v]
            # If no neighbors, diversity is 0
            if not nbrs:
                result[v] = 0
                continue

            # Map old node ids to new (for induced subgraph)
            mapping = {node: idx for idx, node in enumerate(nbrs)}
            # Build induced subgraph: adj list for neighbors only
            sub_n = len(nbrs)
            sub_adj = [[] for _ in range(sub_n)]
            for idx, orig_node in enumerate(nbrs):
                for u in G.adj_list[orig_node]:
                    if u in mapping:
                        sub_adj[idx].append(mapping[u])

            # k-core peeling on induced subgraph
            deg = [len(adj) for adj in sub_adj]
            alive = [True] * sub_n
            queue = deque()
            for i in range(sub_n):
                if deg[i] < k:
                    queue.append(i)
                    alive[i] = False

            # Remove nodes with degree < k
            while queue:
                u = queue.popleft()
                for w in sub_adj[u]:
                    if alive[w]:
                        deg[w] -= 1
                        if deg[w] < k:
                            alive[w] = False
                            queue.append(w)

            # Find connected components in the remaining subgraph
            vis = [False] * sub_n
            cnt = 0
            for i in range(sub_n):
                if alive[i] and not vis[i]:
                    cnt += 1
                    dq = deque()
                    dq.append(i)
                    vis[i] = True
                    while dq:
                        curr = dq.popleft()
                        for nxt in sub_adj[curr]:
                            if alive[nxt] and not vis[nxt]:
                                vis[nxt] = True
                                dq.append(nxt)
            result[v] = cnt

        return result

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
