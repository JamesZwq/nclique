#!/usr/bin/env python3
# Auto-generated for 5506472

STUDENT_ID = "5506472"
STUDENT_NAME = "Yang Song"

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

        def get_kcore_components(adj, nodes, k):
            """
            Finds the number of k-core components in a subgraph induced by nodes.
            """
            # Calculate initial degrees within the subgraph
            deg = {}
            for u in nodes:
                deg[u] = len([v for v in adj.get(u, []) if v in nodes])

            # Initialize queue with nodes having degree less than k
            q = deque([u for u in nodes if deg[u] < k])
            removed = set(q)

            # Iteratively remove nodes and update degrees
            while q:
                u = q.popleft()
                for v in adj.get(u, []):
                    if v in nodes and v not in removed:
                        deg[v] -= 1
                        if deg[v] < k:
                            removed.add(v)
                            q.append(v)

            # Identify k-core nodes
            k_core_nodes = [u for u in nodes if u not in removed]

            # Find connected components within the k-core nodes
            visited = set()
            kcore_components_count = 0

            for u in k_core_nodes:
                if u not in visited:
                    kcore_components_count += 1
                    queue = deque([u])
                    visited.add(u)
                    while queue:
                        v = queue.popleft()
                        for nei in adj.get(v, []):
                            if nei in k_core_nodes and nei not in visited:
                                visited.add(nei)
                                queue.append(nei)

            return kcore_components_count


        n = G.vertex_num
        tau_k = [0] * n

        # Build adjacency list for the entire graph for easier subgraph creation
        graph_adj = G.adj_list

        for v in range(n):
            neighbors = graph_adj[v]
            if not neighbors:
                continue

            # Create the subgraph induced by the neighbors of v
            nbr_nodes = set(neighbors)
            sub_adj = {}
            for u in nbr_nodes:
                sub_adj[u] = [w for w in graph_adj[u] if w in nbr_nodes]


            tau_k[v] = get_kcore_components(sub_adj, list(nbr_nodes), k)

        return tau_k

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
