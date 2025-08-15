#!/usr/bin/env python3
# Auto-generated for 5359871

STUDENT_ID = "5359871"
STUDENT_NAME = "Yuchi Zhang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # Function to get a list of all the connected components formed by
    # a set of nodes (nbs) in a graph G
    @staticmethod
    def getComps(G, nbs):
        comps = []
        visited = set()
        for v in nbs:
            if v not in visited:
                visited.add(v)
                comp = set([v])
                q = deque()
                q.append(v)
                while (q):
                    node = q.popleft()
                    for nb in G.adj_list[node]:
                        if nb not in visited and nb in nbs:
                            q.append(nb)
                            visited.add(nb)
                            comp.add(nb)
                comps.append(comp)
        return comps

    # Function to find the k-core of a subgraph of graph G
    @staticmethod
    def kcore(G, subgraph, k):
        remove = set()
        found = False
        while not found:
            # Remove the nodes with degree less than k from previous iteration
            subgraph = subgraph - remove
            found = True
            # Compute the degree of the node in the subgraph
            for node in subgraph:
                deg = 0
                for nb in G.adj_list[node]:
                    if nb in subgraph:
                        deg += 1
                # If a node has degree less than k, remove it in the next iteration
                if deg < k:
                    remove.add(node)
                    # If a node has been removed in the current iteration,
                    # it means the k-core has not been found yet and we need another iteration
                    found = False

        return subgraph

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
        tau_list = []

        # Compute for all nodes in the graph
        for nbs in G.adj_list:
            # If a node has no neighbours, i.e. is isolated, then it has no
            # neighbourhood induced subgraph and hence tau is 0
            if len(nbs) == 0:
                tau_list.append(0)
                continue
            # Find the k-cores of the neighbourhood induced subgraph given its connected components
            kk_nbs = kCoreBaseStructuralDiversity.kcore(G, set(nbs), k)

            # Find the connected components of the k-cores
            final_kk = kCoreBaseStructuralDiversity.getComps(G, kk_nbs)

            # tau = the number of connected components of the k-cores
            tau = len(final_kk)

            tau_list.append(tau)

        return tau_list


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
