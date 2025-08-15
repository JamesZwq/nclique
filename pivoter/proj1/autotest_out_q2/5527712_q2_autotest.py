#!/usr/bin/env python3
# Auto-generated for 5527712

STUDENT_ID = "5527712"
STUDENT_NAME = "Yuhao Bai"

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
        result = [0] * n
        adj = G.adj_list    # adjacency list provided by the template

        for v in range(n):
            # get the set of v's neighbors
            neighbors = set(adj[v])
            if not neighbors:
                # no neighbors => τ_k(v) = 0
                continue

            # build the neighbor-induced subgraph
            subgraph = {
                u: [w for w in adj[u] if w in neighbors]
                for u in neighbors
            }

            # peel off nodes with degree < k
            core_nodes = kCoreBaseStructuralDiversity._peel_kcore(subgraph, k)
            if not core_nodes:
                # no k-core remains
                continue

            # count connected components among the k-core nodes
            result[v] = kCoreBaseStructuralDiversity._count_components(subgraph, core_nodes)

        return result

    @staticmethod
    def _peel_kcore(adj, k):
        """
        Iteratively remove all vertices of degree < k from adj.
        Return the set of surviving vertices.
        """
        degree = {u: len(adj[u]) for u in adj}
        queue = deque([u for u, d in degree.items() if d < k])
        removed = set(queue)

        while queue:
            u = queue.popleft()
            for w in adj[u]:
                if w not in removed:
                    degree[w] -= 1
                    if degree[w] < k:
                        removed.add(w)
                        queue.append(w)

        return set(adj.keys()) - removed

    @staticmethod
    def _count_components(adj, nodes):
        """
        Count how many connected components exist in adj restricted to 'nodes'.
        """
        seen = set()
        count = 0

        for u in nodes:
            if u not in seen:
                count += 1
                dq = deque([u])
                seen.add(u)
                while dq:
                    x = dq.popleft()
                    for y in adj[x]:
                        if y in nodes and y not in seen:
                            seen.add(y)
                            dq.append(y)

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
