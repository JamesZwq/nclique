#!/usr/bin/env python3
# Auto-generated for 5470173

STUDENT_ID = "5470173"
STUDENT_NAME = "Bin Luo"

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
        sd = [0] * n

        for v in range(n):
            sd[v] = kCoreBaseStructuralDiversity.compute_neighbor_k_core_size(G, v, k)
        return sd


    @staticmethod
    def extract_neighbors(G, v):
        neighbors = G.adj_list[v] # to get all neighbors for each nodes from 0 to n
        if len(neighbors) == 0: # if this node does not have any neighbors, return None
            return None

        neighbor_set = set(neighbors) # to convert to set for faster lookup, time complexity O(d)
        return neighbors, neighbor_set


    @staticmethod
    def discover_neighbor_edges(G, neighbors, neighbor_set):
        neighbor_connections = [] # to save the connection of each nodes, except current node

        # to get each neighbor's edge
        for neighbor in neighbors:
            for adj_of_neighbor in G.adj_list[neighbor]:
                # we don't have to check smaller nodes
                if adj_of_neighbor > neighbor and adj_of_neighbor in neighbor_set: # we need to save their common neighbors
                    neighbor_connections.append((neighbor, adj_of_neighbor))
        # print(neighbor_connections)
        return neighbor_connections


    @staticmethod
    def compute_neighbor_k_core_size(G, v, k):

        # Step1, to get all neighbors for each nodes from 0 to n
        neighbor_result = kCoreBaseStructuralDiversity.extract_neighbors(G, v)
        if neighbor_result is None: # if this node does not have any neighbors, return 0
            return 0
        neighbors, neighbor_set = neighbor_result  # neighbors data

        # Step2, to find edges between neighbors
        neighbor_connections = kCoreBaseStructuralDiversity.discover_neighbor_edges(G, neighbors, neighbor_set)
        # print(neighbor_connections)

        # Step3, to build adjacency list for neighbor-induced subgraph
        adj_list = kCoreBaseStructuralDiversity.build_adjacency_list_from_edges(neighbor_connections)
        # print(adj_list)

        # Step4, to find k-core in neighbor-induced subgraph
        k_core_vertices = kCoreBaseStructuralDiversity.find_k_core_from_adj_list(adj_list, k)
        # print(k_core_vertices)

        # Step5, to calculate connected components using disjoint set
        sd = kCoreBaseStructuralDiversity.count_connected_components_disjoint_set(adj_list, k_core_vertices)
        # print(sd)

        return sd




    @staticmethod
    def build_adjacency_list_from_edges(neighbor_connections):
        if not neighbor_connections: # if no any neighbor connections
            return {}
        adj_list = {}
        # for each edge pairs, add it as key into dictionary
        for u, v in neighbor_connections:
            if u not in adj_list:
                adj_list[u] = []
            if v not in adj_list:
                adj_list[v] = []
        # add value into dict[key]
        for u, v in neighbor_connections:
            adj_list[u].append(v)
            adj_list[v].append(u)
        return adj_list



    @staticmethod
    def find_k_core_from_adj_list(adj_list, k):
        if not adj_list: # if adjacent list is all empty
            return []

        degrees = {}
        for vertex in adj_list:
            degrees[vertex] = len(adj_list[vertex])

        queue = []
        removed = set()
        for vertex_id in degrees:
            vertex_degree = degrees[vertex_id]
            if vertex_degree < k:
                queue.append(vertex_id)
                removed.add(vertex_id)

        # as lecture said, we can delete node iteratively
        queue = deque(queue) # change it to deque
        while queue:
            vertex = queue.popleft() # O(1)

            for neighbor in adj_list[vertex]:
                if neighbor not in removed:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        queue.append(neighbor) # O(1)
                        removed.add(neighbor)

        k_core_vertices = list(set(adj_list.keys()) - removed) # return nodes in the k_core
        return k_core_vertices




    # As lecutre said, use disjoint set to find number of connected subgraphs
    @staticmethod
    def count_connected_components_disjoint_set(adj_list, vertices):
        if not vertices:
            return 0

        # to initialize disjoint set data structure
        parent = {}
        setsize = {}
        for e in vertices:
            parent[e] = e # to set each vertex as its own parent initially
            setsize[e] = 1 # to set size of each component to 1

        vertices_set = set(vertices) # to convert to set for faster search

        def Find(e):
            # to find root with simple path following (like C version)
            while parent[e] != e:
                e = parent[e]
            return e

        def Union(i, j):
            # to union by size exactly like the C implementation
            i = Find(i)
            j = Find(j)
            if i == j:
                return False # already connected

            if setsize[i] < setsize[j]:
                # to link smaller tree to larger one
                setsize[j] += setsize[i]
                parent[i] = j
            else:
                # to link smaller tree to larger one
                setsize[i] += setsize[j]
                parent[j] = i

            return True

        component_count = len(vertices)  # to count total components initially

        # to process edges and merge components
        for vertex in vertices:
            if vertex in adj_list:
                for neighbor in adj_list[vertex]:
                    if neighbor in vertices_set and vertex < neighbor:
                        if Union(vertex, neighbor): # to merge if different components
                            component_count -= 1 # to decrease component count

        return component_count

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
