#!/usr/bin/env python3
# Auto-generated for 5571439

STUDENT_ID = "5571439"
STUDENT_NAME = "Zebu Chen"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

from collections import deque
from collections import defaultdict
import heapq
# Time complexity:O(V+E+V*d^2)
# V and E is the total number of vertices and edges
# d is the maximum degree in the graph
# At worst case like complete graph,the time complexity will be O(V^3)
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
        # # TODO
        n = G.vertex_num

        sd = [0] * n
        if k <= 0 or n == 0:
            return sd
        core_number_list=kCoreBaseStructuralDiversity.compute_core_number(G)
        # print(core_number_list)
        # all nodes whose core number ≥ k+1
        # on the original graph


        for i in range(n):
            # skip those vertices whose degree < k
            if core_number_list[i]<k:
                continue
            # construct a new subgraph of all vertices whose core number>k.
            subgraph_nodes = set()
            for nb in G.adj_list[i]:
                if core_number_list[nb] >= k:
                    subgraph_nodes.add(nb)
            if not subgraph_nodes:
                continue
            # k:vertex , v:the adj_list of the vertex
            new_adjacency_list = defaultdict(list)
            for i_qual_nb in subgraph_nodes:
                for both_nb in G.adj_list[i_qual_nb]:
                    if both_nb in subgraph_nodes:
                        new_adjacency_list[i_qual_nb].append(both_nb)

            # prune to get a new subgraph of only vertices with a core number≥k
            n_subgraph=kCoreBaseStructuralDiversity.k_core_subgraph_prune(new_adjacency_list, k)
        # Sort nodes by degree in descending order
            sd[i]=(kCoreBaseStructuralDiversity.calculate_connected_components(n_subgraph))
        # print(sd)
        return sd


    @staticmethod
    # input : the subgraph of all neighbors of the chosen vertex whose core number ≥ k+1 , and core number k
    # return : the new subgraph with all vertices whose core number ≥ k
    def k_core_subgraph_prune(subgraph,k):
        # make a copy,then use set to save more time
        cur = {node: set(neighs) for node, neighs in subgraph.items()}
        q=deque()
        for node in list(cur.keys()):
            if len(cur[node])<k:
                q.append(node)
        while q:
            remove_node=q.popleft()
            if remove_node not in cur:
                continue
            for remove_nb in cur[remove_node]:
                if remove_nb in cur:
                # cur[remove_nb].remove(remove_node)
                    cur[remove_nb].discard(remove_node)
                if len(cur[remove_nb])<k and remove_nb not in q:
                    q.append(remove_nb)
            del cur[remove_node]
        return cur


    @staticmethod
    # BFS to calculate the number of the connected components
    # input: default dict,k is vertex ,v is its adjacent list
    def calculate_connected_components(subgraph_dict):
        visited = set()
        connected_components_num = 0

        for key in subgraph_dict:
            if key not in visited:
                q =deque()
                q.append(key)
                visited.add(key)
                while q:
                    cur_node=q.popleft()
                    for nb in subgraph_dict[cur_node]:
                        if nb not in visited:
                            visited.add(nb)
                            q.append(nb)
                connected_components_num += 1
        # return the number of connected components
        return connected_components_num



    @staticmethod
    def compute_core_number(G):
        n = G.vertex_num
        degree_list = [len(G.adj_list[i]) for i in range(n)]
        core_list = degree_list[:]
        # k: the value of degree
        current_degree = defaultdict(list)
        max_degree = max(degree_list)

        # put all vertices in the position of its degree
        for v in range(n):
            current_degree[degree_list[v]].append(v)

        # ascent order of degree
        # 'remove' nodes from the original graph
        for d in range(max_degree + 1):
            # vertices of same degree in the same batch
            while current_degree[d]:
                v = current_degree[d].pop()
                for u in G.adj_list[v]:
                    # only if the degree of its neighbor ≥ the vertex
                    # can it means that the neighbor hasn't been proceeded
                    if degree_list[u] > d:
                        # 'remove' v from u's neighbors
                        degree_list[u] -= 1
                        current_degree[degree_list[u]].append(u)
                # mark the core number as the current degree
                core_list[v] = d
        return core_list



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
