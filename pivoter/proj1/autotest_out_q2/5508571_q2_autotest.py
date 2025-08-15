#!/usr/bin/env python3
# Auto-generated for 5508571

STUDENT_ID = "5508571"
STUDENT_NAME = "Jiahao Zhang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
import collections
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
        sd = [0] * n  # list to store τ_k(v) for all vertices

        for v in range(n):
            # get neighbors of v
            neighbors_of_v = G.adj_list[v]

            # if no neighbors, τ_k(v) = 0
            if not neighbors_of_v:
                sd[v] = 0
                continue

            # build the neighbour-induced subgraph H
            neighbors_set = set(neighbors_of_v)
            nbr_adj = {}
            # for each neighbor of v, create an empty list of adjacent neighbors
            for node in neighbors_of_v:
                nbr_adj[node] = []

            for u1 in neighbors_of_v:
                for u2 in G.adj_list[u1]:
                    if u1 < u2 and u2 in neighbors_set:
                        nbr_adj[u1].append(u2)
                        nbr_adj[u2].append(u1)

            # map each neighbour to a temporary ID
            sorted_neighbors = sorted(neighbors_of_v)
            # initialize the mapping dictionary
            original_to_temp = {}
            # assign a temporary ID to each neighbor
            i = 0
            for node in sorted_neighbors:
                original_to_temp[node] = i
                i += 1
            # build the reverse mapping from temporary ID back to original node
            temp_to_original = {}
            for node, i in original_to_temp.items():
                temp_to_original[i] = node

            # ompute core numbers in H
            nbr_core_nums = kCoreBaseStructuralDiversity._cal_core_nums_for_subgraph(nbr_adj, original_to_temp, temp_to_original)

            # select k-core candidate vertices
            candidates = []
            for node, cnum in nbr_core_nums.items():
                if cnum >= k:
                    candidates.append(node)
            if not candidates:
                sd[v] = 0
                continue

            # build k-core-induced subgraph of H
            hk_adj = {}
            for node in candidates:
                hk_adj[node] = []
            for u1 in candidates:
                for u2 in nbr_adj[u1]:
                    if u2 in hk_adj:
                        hk_adj[u1].append(u2)
            # remove duplicate edges
            for node in hk_adj:
                # get the current neighbor list
                neighbors = hk_adj[node]

                # remove duplicates by converting to a set
                unique_neighbors = set(neighbors)

                # convert the set back to a list
                unique_list = list(unique_neighbors)

                # assign the deduplicated list back
                hk_adj[node] = unique_list

            # count connected components in the k-core subgraph
            visited = set()
            num_components = 0
            for start in candidates:
                if start not in visited:
                    num_components += 1
                    queue = collections.deque([start])
                    visited.add(start)
                    while queue:
                        u = queue.popleft()
                        for w in hk_adj[u]:
                            if w not in visited:
                                visited.add(w)
                                queue.append(w)

            sd[v] = num_components

        return sd

    @staticmethod
    def _cal_core_nums_for_subgraph(adj_list, orig_to_temp, temp_to_orig):
        """
        Compute core numbers for a subgraph via a bucket-based peeling algorithm.
        """
        n = len(adj_list)
        if n == 0:
            return {}

        #convert to temporary-ID adjacency list and degrees
        temp_adj = [[] for _ in range(n)]
        degrees = [0] * n
        for orig_u, nbrs in adj_list.items():
            tu = orig_to_temp[orig_u]
            for orig_v in nbrs:
                tv = orig_to_temp[orig_v]
                temp_adj[tu].append(tv)
            degrees[tu] = len(temp_adj[tu])

        core = [0] * n
        max_deg = max(degrees) if degrees else 0

        # bucket sort vertices by degree
        V = [0] * n
        pos = [0] * n
        bin_counts = [0] * (max_deg + 1)

        # count how many vertices have each degree
        for d in degrees:
            bin_counts[d] += 1

        # 'start' will track the beginning index for each degree bucket in V
        start = 0

        #convert bin_counts from counts to starting positions
        for d in range(max_deg + 1):
            count = bin_counts[d]
            bin_counts[d] = start
            start += count

        # plaace each vertex into V at the index corresponding to its degree
        for v_temp in range(n):
            d = degrees[v_temp]
            # record the position of v_temp in V
            pos[v_temp] = bin_counts[d]
            V[bin_counts[d]] = v_temp
            # move the bucket start so next same‑degree vertex goes into the next slot
            bin_counts[d] += 1

        # reset bin_counts so that bin_counts[d] again points to the first index of degree‑d bucket
        for d in range(max_deg, 0, -1):
            bin_counts[d] = bin_counts[d - 1]
        bin_counts[0] = 0
        
        i = 0
        while i < n:
            # get the temporary vertex ID at position i
            v_temp = V[i]
            
            # record its current degree as its core number
            deg_v_temp = degrees[v_temp]
            core[v_temp] = deg_v_temp

            # get the list of neighbors (temporary IDs)
            neighbor_list = temp_adj[v_temp]
            
            # iterate through each neighbor
            j = 0
            while j < len(neighbor_list):
                nbr = neighbor_list[j]
                
                # check if this neighbor has higher degree
                deg_nbr = degrees[nbr]
                if deg_nbr > deg_v_temp:
                    # du = neighbor's degree
                    du = deg_nbr
                    # pu = neighbor's position in V
                    pu = pos[nbr]
                    # pw = start index of the bucket for degree du
                    pw = bin_counts[du]

                    # swap V[pu] and V[pw]
                    temp_vertex = V[pu]
                    V[pu] = V[pw]
                    V[pw] = temp_vertex

                    # update positions in pos[]
                    pos[V[pu]] = pu
                    pos[V[pw]] = pw

                    # advance the bucket start for degree du
                    bin_counts[du] = bin_counts[du] + 1

                    # decrement the stored degree of the neighbor
                    degrees[nbr] = degrees[nbr] - 1

                j += 1

            i += 1
        # build the final mapping from original IDs to core numbers
        result = {}
        for i in range(n):
            orig_id = temp_to_orig[i]
            core_num = core[i]
            result[orig_id] = core_num

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
