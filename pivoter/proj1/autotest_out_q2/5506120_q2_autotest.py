#!/usr/bin/env python3
# Auto-generated for 5506120

STUDENT_ID = "5506120"
STUDENT_NAME = "Sho Tazume"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
import heapq    # import for heap version peeling algorithm
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
        # TODO
        # get the numver of vertices
        n = G.vertex_num
        # initialise tau list
        tau = [0] * n

        # for each vertex v do the step below
        #    step 1. get the list of neighbours nbr(v)
        #    step 2. construct the induced subgraph of nbr(v)
        #    step 3. compute the k-core of this subgraph
        #    step 4. count the number of connected commponents in the k-core

        # loop all vertices one by one
        for v in range(n):
            # step 1
            # get neighbours of vertex v
            nbrs = G.adj_list[v]    # nbr(v)
            
            # construct the induced subgraph of nbr(v)
            nodes_in_subgraph = set(nbrs)   # that is neighbours
            edges_in_subgraph = []          # initially no edges

            # construct the neighbour induced subgraph
            # pick one neighbour of v one by one
            for u in nodes_in_subgraph:
                # pick the neighbours of u one by one
                for w in G.adj_list[u]:
                    # if w is a neighbour of v, add the edge
                    if w in nodes_in_subgraph:
                        # to avoid duplicate edges, add only if u < w
                        if u < w:
                            edges_in_subgraph.append((u, w))

            # step 2
            # create the induced subgraph
            induced_adj = {u: [] for u in nodes_in_subgraph}
            for u, w in edges_in_subgraph:
                induced_adj[u].append(w)
                induced_adj[w].append(u)

            # step 3
            # extract k-core
            # k_core_nodes = _extract_k_core(induced_adj, k)
            k_core_nodes = kCoreBaseStructuralDiversity._extract_k_core(induced_adj, k)
            
            # step 4
            # count the number of connected components in the k-core
            # tau[v] = _count_num_of_connected_components(induced_adj, k_core_nodes)
            tau[v] = kCoreBaseStructuralDiversity._count_num_of_connected_components(induced_adj, k_core_nodes)
            
        return tau
            
                            
    # step 3
    # calculate k core by using decomposition
    # use priority queue to extract k-core because it's more efficient
    @staticmethod
    def _extract_k_core(adj_list: dict, k: int)-> set:
        """
        extract k-core from the adjacency list of the graph with Global-view peeling algorithm
        """
        # degree dictionary: key is node, value is degree
        degree = {u: len(adj_list[u]) for u in adj_list}
        # heap containing node and its degree
        heap = [(deg, u) for u, deg in degree.items()]
        # convert to heap
        heapq.heapify(heap)     # the degree with the smallest degree comes first
        # set to record visited nodes
        visited = set()

        while heap:
            # extract the node with the smallest degree
            deg_u, u = heapq.heappop(heap)
            # if we already checked the node, skip it
            if u in visited:
                continue
            # if the smallest remaining degree is already >= k, remaining nodes are part of k-core
            if deg_u >= k:
                break
            # mark the node as visited
            visited.add(u)

            # decrease the degree of the neighbours of u with higher degree than u by 1
            for v in adj_list[u]:
                # if v is not visited, and v's degree is higher than u decrease its degree by1
                # if v not in visited:
                if v not in visited and degree[v] > degree[u]:
                    degree[v] -= 1
                    # update the heap with the new degree of v
                    heapq.heappush(heap, (degree[v], v))

        k_core = set()
        # return the set of nodes whose degree is more than or equal to k
        for u in adj_list:
            # check the vertex is not visited
            if u not in visited:
                # check the degree of u is more than or equal to k
                if degree[u] >= k:
                    k_core.add(u)
        return k_core


    # step 4
    @staticmethod
    def _count_num_of_connected_components(adj_list: dict, nodes: set) -> int:
        parent = {u: u for u in nodes}

        def find(u: int) -> int:
            """
            find the root of node u with path compression
            """
            # go up the tree until we find the root
            while parent[u] != u:
                # path compression
                parent[u] = parent[parent[u]]  
                u = parent[u]
            return u

        def union(u, v):
            """
            union the two components nodes u and v
            """
            # find the roots of u and v
            pu, pv = find(u), find(v)
            if pu != pv:
                parent[pu] = pv

        # apply union operations for each edge between nodes
        for u in nodes:
            for v in adj_list[u]:
                if v in nodes:
                    union(u, v)

        # count the number of distinct roots
        roots = set(find(u) for u in nodes)
        return len(roots)

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
