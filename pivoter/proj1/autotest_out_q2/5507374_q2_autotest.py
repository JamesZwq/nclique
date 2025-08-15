#!/usr/bin/env python3
# Auto-generated for 5507374

STUDENT_ID = "5507374"
STUDENT_NAME = "(Albert) Simeng Li"

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

        # De-duplicate all neighbor lists
        for v in range(n):
            G.adj_list[v] = list(set(G.adj_list[v]))

        find_nbrs = lambda v: G.adj_list[v]
        sd = [0] * n

        # For each vertex v, treat its neighboring subgraphs individually
        for v in range(n):
            Ns = find_nbrs(v)
            if not Ns:
                sd[v] = 0
                continue

            Nset = set(Ns)
            # Calculate the degree of each u in the subgraph
            degrees = {}
            for u in Ns:
                # Ensure that neighbors are not duplicated with set(find_nbrs(u))
                cnt = 0
                for i in set(find_nbrs(u)):
                    if i in Nset:
                        cnt += 1
                degrees[u] = cnt

            # k-core pruning on subgraphs
            inG = {u: True for u in Ns}
            dq = deque(u for u, d in degrees.items() if d < k)
            while dq:
                u = dq.popleft()
                if not inG[u]:
                    continue
                inG[u] = False
                for v in set(find_nbrs(u)):
                    if v in Nset and inG[v]:
                        degrees[v] -= 1
                        if degrees[v] < k:
                            dq.append(v)

            # Count the number of connected components of the remaining nodes
            seen = set()
            count = 0
            for u in Ns:
                if inG.get(u, False) and u not in seen:
                    count += 1
                    sq = deque([u])
                    seen.add(u)
                    while sq:
                        x = sq.popleft()
                        for i in set(find_nbrs(x)):
                            if i in Nset and inG.get(i, False) and i not in seen:
                                seen.add(i)
                                sq.append(i)

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
