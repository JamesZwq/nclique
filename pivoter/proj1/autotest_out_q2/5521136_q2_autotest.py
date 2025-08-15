#!/usr/bin/env python3
# Auto-generated for 5521136

STUDENT_ID = "5521136"
STUDENT_NAME = "Xi Yan"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################


class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Main processing function to compute τ_k(v) for each node:
        For each node v, extract the induced subgraph formed by its neighbors,
        retain only the k-core portion, and count the number of connected components
        as τ_k(v).
        
        Parameters:
        - G: The input graph object that contains adjacency list and other structural information.
        - k: The core decomposition threshold (minimum degree constraint).

        Returns a list of size G.vertex_num, 
        where result[v] represents the structural diversity τ_k(v)
        """
        n = G.vertex_num
        result = [0] * n   # Initialize result list

        # Convert adjacency lists to sets for efficient intersection operations
        neigh_sets = [set(nei) for nei in G.adj_list]

        # Iterate over each node v to compute its τ_k(v)
        for v, nbrs in enumerate(neigh_sets):
            # If the number of available neighbors is less than k, leave 0 and skip directly
            if len(nbrs) < k:
                continue

            # 1) Compute the degree of each node in the subgraph: number of common neighbors between u and v
            deg_map = {}
            for u in nbrs:
                deg_map[u] = len(neigh_sets[u] & nbrs)

            # 2) Initialize peeling queue
            # collect all nodes with degree less than k
            peel = deque([u for u, d in deg_map.items() if d < k])
            deleted = set(peel)

            # 3) Iterative peeling
            # for each removed node, update the degrees of its remaining neighbors
            while peel:
                u = peel.popleft()
                # Only iterate over neighbors of u that are in the core subgraph
                for w in neigh_sets[u] & nbrs:
                    if w in deleted:
                        continue
                    deg_map[w] -= 1
                    # Once the degree drops below k, add the node to the deleted set and continue peeling
                    if deg_map[w] < k:
                        deleted.add(w)
                        peel.append(w)

            # 4) The remaining nodes form the k-core
            # count its connected components
            core_nodes = nbrs - deleted
            result[v] = kCoreBaseStructuralDiversity._count_components(G, core_nodes)

        return result

    @staticmethod
    def _count_components(G, nodes):
        """
        Helper function: Count the number of connected components in the subgraph
        induced by the given set of nodes.

        Parameters:
        - G: the original graph object
        - nodes: a set of nodes forming the induced subgraph

        Returns:
        - comps: number of connected components (integer)

        Count the number of connected components in the subgraph 
        defined by the given set of nodes.
        """
        # No components in an empty graph
        if not nodes:
            return 0

        seen = set()  # Track visited nodes
        comps = 0  # Component counter
        
        # DFS to find each connected component
        for u in nodes:
            if u in seen:
                continue
            # Use a list as a stack to perform depth-first traversal
            stack = [u]
            seen.add(u)
            while stack:
                cur = stack.pop()
                for w in G.adj_list[cur]:
                    if w in nodes and w not in seen:
                        seen.add(w)
                        stack.append(w)
            # Each DFS marks one connected component
            comps += 1

        return comps

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
