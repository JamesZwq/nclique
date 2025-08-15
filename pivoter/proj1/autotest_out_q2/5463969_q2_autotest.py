#!/usr/bin/env python3
# Auto-generated for 5463969

STUDENT_ID = "5463969"
STUDENT_NAME = "Yicheng Xu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def build_induced_subgraph(G, neighbors):
        """
        Constructs the induced subgraph of a given node's neighborhood.
        Only includes edges between neighbors of the node.
        """
        sub_adj = defaultdict(set)
        for u in neighbors:
            for v in G.adj_list[u]:
                if v in neighbors:
                    sub_adj[u].add(v)
                    sub_adj[v].add(u)
        return sub_adj

    @staticmethod
    def k_core_peeling(sub_adj, k):
        """
        Performs k-core peeling: iteratively removes nodes with degree < k.
        Returns the remaining node set that satisfies the k-core condition.
        """
        degree = {node: len(neighs) for node, neighs in sub_adj.items()}
        q = deque([node for node, deg in degree.items() if deg < k])
        while q:
            node = q.popleft()
            for neigh in sub_adj[node]:
                if degree[neigh] >= k:
                    degree[neigh] -= 1
                    if degree[neigh] == k - 1:
                        q.append(neigh)
            degree[node] = -1
        return {node for node, deg in degree.items() if deg >= k}

    @staticmethod
    def count_components(sub_adj, valid_nodes):
        """
        Counts the number of connected components in a given valid subgraph.
        Uses standard DFS for component labeling.
        """
        visited = set()
        count = 0
        for node in valid_nodes:
            if node not in visited:
                count += 1
                stack = [node]
                visited.add(node)
                while stack:
                    u = stack.pop()
                    for v in sub_adj[u]:
                        if v in valid_nodes and v not in visited:
                            visited.add(v)
                            stack.append(v)
        return count

    @staticmethod
    def process(G, k):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int (k-core level)

        Returns
        -------
        List[int]  # τ_k(v) for all v

        Notes
        -----
        - For each vertex v, we construct the neighbor-induced subgraph G[N(v)]
        - Compute the number of connected components that are valid k-cores
        - A valid k-core is a connected subgraph where each node has degree >= k
        - This function handles edge cases like isolated nodes and k=0
        - Also handles globally edgeless graphs (no neighbors for any vertex)

        Time Complexity
        ---------------
        Let d = average degree of a node
        For each node v:
        - Building induced subgraph: O(d^2)
        - Peeling process (k-core): O(d^2)
        - Component counting: O(d^2)
        Overall: O(n * d^2) ~ linear for sparse graphs
        """
        n = G.vertex_num
        tau = [0] * n

        # Fast check for globally edgeless graph
        if all(len(G.adj_list[v]) == 0 for v in range(n)):
            return tau  # all tau[v] = 0

        for v in range(n):
            neighbors = set(G.adj_list[v])

            # Handle isolated node
            if not neighbors:
                tau[v] = 0
                continue

            # Build the neighbor-induced subgraph
            sub_adj = kCoreBaseStructuralDiversity.build_induced_subgraph(G, neighbors)

            # Handle k = 0: count all components in the neighbor-induced subgraph
            if k == 0:
                tau[v] = kCoreBaseStructuralDiversity.count_components(sub_adj, neighbors)
            else:
                # Peeling process: remove nodes with degree < k
                core_nodes = kCoreBaseStructuralDiversity.k_core_peeling(sub_adj, k)
                if core_nodes:
                    tau[v] = kCoreBaseStructuralDiversity.count_components(sub_adj, core_nodes)
                else:
                    tau[v] = 0  # No valid k-core found

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
