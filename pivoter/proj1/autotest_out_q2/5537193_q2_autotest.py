#!/usr/bin/env python3
# Auto-generated for 5537193

STUDENT_ID = "5537193"
STUDENT_NAME = "Hongpeng Hua"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################
from collections import deque

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
        #  Compatible with different adjacency list fields
        if hasattr(G, "adj_list"):
            adjacency_list = G.adj_list
        elif hasattr(G, "adj_list_out"):
            adjacency_list = G.adj_list_out
        else:
            raise AttributeError("Graph object has no adj_list attribute") # If neither of the two attributes is present, an error will be reported。
        n = G.vertex_num
        tau = [0] * n   # Prepare a list of length n to store the result of each vertex

        # Calculate τ_k(v) for each vertex v
        for v in range(n):
            Neighbours = adjacency_list[v]
            deg_v = len(Neighbours)
            if deg_v == 0:
                tau[v] = 0
                continue

            # Neighbor set marker
            belong = [False] * n
            for u in Neighbours:
                belong[u] = True

            # Count the degrees in the neighbour induction subgraph
            subdeg = {}
            for u in Neighbours:
                cnt = 0
                for w in adjacency_list[u]:
                    if belong[w]:     # Only count the edges between neighbours
                        cnt += 1
                subdeg[u] = cnt

            if k == 0:  # When k=0, there is no need to peel; directly count the connected components
                alive = {u: True for u in Neighbours}
            else:
                # k-core
                alive = {u: True for u in Neighbours}
                q = deque([u for u in Neighbours if subdeg[u] < k])
                while q:
                    u = q.popleft()
                    if not alive[u]:
                        continue
                    alive[u] = False
                    # Update the degree of each of its neighbours
                    for w in adjacency_list[u]:
                        if belong[w] and alive.get(w, False):
                            subdeg[w] -= 1
                            if subdeg[w] < k:
                                q.append(w)

            #  Count the number of connected components of the remaining "alive" points
            visited = set()
            component_count = 0
            for u in Neighbours:
                if alive.get(u, False) and u not in visited:
                    component_count += 1
                    dq = deque([u])
                    visited.add(u)   # BFS marks all points within this block to prevent double counting
                    while dq:
                        x = dq.popleft()
                        for y in adjacency_list[x]:
                            if belong[y] and alive.get(y, False) and y not in visited:
                                visited.add(y)
                                dq.append(y)

            tau[v] = component_count   # store the calculated number of components in the result list

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
