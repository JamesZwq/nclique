#!/usr/bin/env python3
# Auto-generated for 5520221

STUDENT_ID = "5520221"
STUDENT_NAME = "Yuting Liang"

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
        n = G.vertex_num
        sd = [0] * n

        # iterate over each vertex in the graph
        for v in range(n):
            neighbors = G.adj_list[v]  # The neighbor list of vertex v

            # if v has no neighbors
            if not neighbors:
                sd[v] = 1 if k == 0 else 0   #vertex with no neighbors has 1 0-core or 0 otherwise
                continue

            # build the neighbor-induced subgraph of v
            # map from original vertex id to index in neighbors
            idx_map = {}
            for i, u in enumerate(neighbors):
                idx_map[u] = i

            # initialize adjacency list for the subgraph
            subgraph = [[] for _ in range(len(neighbors))]

            # fill subgraph adjacency list with valid edges between neighbors
            for i, u in enumerate(neighbors):
                for w in G.adj_list[u]:
                    if w in idx_map and w != u:  #  no self-loops
                        subgraph[i].append(idx_map[w])

            # iteratively perform k-core peeling to remove vertices with degree < k
            deg = [len(adj) for adj in subgraph]  # current degree in the subgraph
            active_v = set(range(len(neighbors)))  # set of not yet removed nodes
            queue = deque([i for i in active_v if deg[i] < k])  # initially removable nodes

            while queue:
                curr = queue.popleft()
                if curr not in active_v:
                    continue
                active_v.remove(curr)
                for w in subgraph[curr]:
                    if w in active_v:
                        deg[w] -= 1  # decrease degree by 1
                        if deg[w] < k:
                            queue.append(w)

            # count number of connected components in the remaining subgraph
            visited = set()
            count = 0  # bumber of connected components

            for u in active_v:
                if u in visited:
                    continue
                # DFS from u to find connected component
                stack = [u]
                visited.add(u)
                while stack:
                    curr = stack.pop()
                    for next_v in subgraph[curr]:
                        if next_v in active_v and next_v not in visited:
                            visited.add(next_v)
                            stack.append(next_v)
                count += 1  # component found

            # save result
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
