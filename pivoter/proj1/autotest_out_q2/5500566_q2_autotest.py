#!/usr/bin/env python3
# Auto-generated for 5500566

STUDENT_ID = "5500566"
STUDENT_NAME = "Wei Wang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
import collections
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
        n = G.vertex_num            # total number of vertices in the graph
        sd = [0] * n                # output list: sd[v] = τ_k(v) for each v

        for v in range(n):
            # Step 1: Find all neighbors of vertex v (the neighbour set N(v))
            nbrs = G.adj_list[v]
            if not nbrs:
                # If v has no neighbors, its neighbour-induced subgraph is empty
                sd[v] = 0
                continue

            # Step 2: Build the neighbour-induced subgraph H for v
            nbr_set = set(nbrs)  # for O(1) lookups
            # sub_adj stores the adjacency list of H, i.e. induced by N(v)
            sub_adj = {node: [] for node in nbrs}
            for u in nbrs:
                for w in G.adj_list[u]:
                    # only include edges where both endpoints are in N(v)
                    if w in nbr_set and w > u:
                        sub_adj[u].append(w)
                        sub_adj[w].append(u)

            # Step 3: Map nodes in H to temporary integer IDs for fast bucket peeling
            nodes_sorted = sorted(nbrs)
            orig2temp = {orig: idx for idx, orig in enumerate(nodes_sorted)}  # original id to temp id
            temp2orig = {idx: orig for orig, idx in orig2temp.items()}        # temp id to original id

            # Step 4: Compute core number for each node in H
            core_nums = kCoreBaseStructuralDiversity._calc_core_numbers(sub_adj, orig2temp, temp2orig)

            # Step 5: Select nodes in H whose core number >= k to form the k-core
            k_core_nodes = [node for node, cn in core_nums.items() if cn >= k]
            if not k_core_nodes:
                # No k-core exists in the neighbour-induced subgraph of v
                sd[v] = 0
                continue

            # Step 6: Build adjacency for the induced k-core subgraph of H
            # Only keep nodes and edges between k-core nodes
            kcore_adj = {node: [w for w in sub_adj[node] if w in k_core_nodes] for node in k_core_nodes}

            # Step 7: Count connected components in the k-core subgraph of H using BFS
            comp_count = 0              # number of connected components
            seen = set()                # set of already visited nodes
            for start in k_core_nodes:
                if start not in seen:
                    comp_count += 1
                    q = deque([start])  # use a queue for BFS
                    seen.add(start)
                    while q:
                        cur = q.popleft()
                        for nxt in kcore_adj[cur]:
                            if nxt not in seen:
                                seen.add(nxt)
                                q.append(nxt)
            # Store result for current vertex v
            sd[v] = comp_count

        return sd

    @staticmethod
    def _calc_core_numbers(adj_dict, orig2temp, temp2orig):
        """
        Calculate the core number for every node in a given subgraph using bucket-based peeling algorithm.

        Parameters
        ----------
        adj_dict : dict
            Adjacency list of the subgraph, keys are original node ids
        orig2temp : dict
            Mapping from original node id to temp integer id (for array access)
        temp2orig : dict
            Reverse mapping from temp id to original id

        Returns
        -------
        dict : Mapping from original node id to its core number in the subgraph
        """
        n = len(adj_dict)               # number of nodes in subgraph
        if n == 0:
            return {}

        # Convert to arrays using temp ids for efficiency
        temp_adj = [[] for _ in range(n)]      # temp_id adjacency list
        degree = [0] * n                      # degree for each temp_id
        for orig_u, nbrs in adj_dict.items():
            tu = orig2temp[orig_u]
            for orig_v in nbrs:
                tv = orig2temp[orig_v]
                temp_adj[tu].append(tv)
            degree[tu] = len(temp_adj[tu])

        core = [0] * n                        # core numbers for temp ids
        max_deg = max(degree) if degree else 0

        # Bucket sort: V contains all nodes, pos gives their position in V
        V = [0] * n
        pos = [0] * n
        bin_starts = [0] * (max_deg + 1)      # start position of each degree in V

        # Count how many nodes have each degree
        for d in degree:
            bin_starts[d] += 1
        start = 0
        for d in range(max_deg + 1):
            count = bin_starts[d]
            bin_starts[d] = start
            start += count
        # Place each node into V according to its degree
        for v_temp in range(n):
            d = degree[v_temp]
            pos[v_temp] = bin_starts[d]
            V[bin_starts[d]] = v_temp
            bin_starts[d] += 1
        # Restore bin_starts to point to start of each degree bucket
        for d in range(max_deg, 0, -1):
            bin_starts[d] = bin_starts[d - 1]
        bin_starts[0] = 0

        # Main peeling loop: repeatedly remove lowest-degree nodes
        for i in range(n):
            v_temp = V[i]
            core[v_temp] = degree[v_temp]
            for u_temp in temp_adj[v_temp]:
                if degree[u_temp] > degree[v_temp]:
                    du = degree[u_temp]
                    pu = pos[u_temp]
                    pw = bin_starts[du]
                    # Swap u_temp into the current start of its degree bucket
                    V[pu], V[pw] = V[pw], V[pu]
                    pos[V[pu]], pos[V[pw]] = pu, pw
                    bin_starts[du] += 1
                    degree[u_temp] -= 1

        # Convert results from temp ids back to original node ids
        result = {}
        for idx in range(n):
            orig_id = temp2orig[idx]
            result[orig_id] = core[idx]
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
