#!/usr/bin/env python3
# Auto-generated for 5509272

STUDENT_ID = "5509272"
STUDENT_NAME = "Yifan Tan"

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

        for v in range(n): # O(N)
            neighbors = G.adj_list[v]
            # No neighbors means no k-core in the induced subgraph
            if not neighbors:
                sd[v] = 0
                continue

            # Boolean array to mark nodes currently in the core
            in_core = [False] * n
            for u in neighbors: # O(degv) Worst # O(N), degv is the degree of node v
                in_core[u] = True

            # Compute local degrees within neighbor-induced subgraph
            # O(∑_{u∈N(v)} degu) Worst O(N^2), N(v) denotes the set of neighboring nodes of v
            local_deg = [0] * n
            for u in neighbors: # O(degv)
                cnt = 0
                for w in G.adj_list[u]: # O(degu)
                    if in_core[w]:
                        cnt += 1
                local_deg[u] = cnt

            # k-core peeling: iteratively remove nodes with degree < k
            # O(∑_{u∈N(v)} degu) Worst O(N^2)
            queue = deque([u for u in neighbors if local_deg[u] < k])
            while queue:
                x = queue.popleft()
                if not in_core[x]:
                    continue
                in_core[x] = False
                # Update degrees of neighbors
                for y in G.adj_list[x]:
                    if in_core[y]:
                        local_deg[y] -= 1
                        if local_deg[y] < k:
                            queue.append(y)

            # DFS Count connected components in the remaining k-core
            # O(∑_{u∈N(v)} degu) Worst O(N+M)
            comp_count = 0
            visited = [False] * n
            for u in neighbors: # O(degu)
                if in_core[u] and not visited[u]:
                    comp_count += 1
                    stack = [u]
                    visited[u] = True
                    while stack:
                        node = stack.pop()
                        for w in G.adj_list[node]:
                            if in_core[w] and not visited[w]:
                                visited[w] = True
                                stack.append(w)

            sd[v] = comp_count

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
