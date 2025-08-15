#!/usr/bin/env python3
# Auto-generated for 5573658

STUDENT_ID = "5573658"
STUDENT_NAME = "Taliesen Barker"

# ======= 学生代码 =======
"""
My solution to this question consists of a sequence of calls to optimal implementations of well-known graph algorithms.
For each node, we begin by computing the neighbour-induced subgraph. With the subgraph we then calculate the core number
for each node - to do this we use bin sort to sort nodes by degree, and a flat array implementation to do core peeling.
Finally, we reconstruct the connected components based on the subgraph and the core numbers assigned in the previous 
step, using a Discrete Set data structure.
In total, the worst case time complexity for this solution is O(n^3 + n*m*log(n)), though on average it should be much better
than this.
Note: on an old laptop with the sample tests my code takes up to 5 seconds to run."""

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
        return kCoreBaseStructuralDiversity.kcore_structural_diversity(G, k)

    @staticmethod
    def kcore_structural_diversity(graph, k):
        """determines the kcore structural diversity for
        each node in the given graph"""
        # total O(n*(n^2+m+n+n+m*log(n)))->O(n^3+nm*log(n))
        n = graph.vertex_num
        edge_list = graph.adj_list
        result = [0] * n
        for node in range(n):  # O(n)
            subgraph = kCoreBaseStructuralDiversity.get_neighbour_subgraph(edge_list, node)  # O(n^2)
            core_numbers = kCoreBaseStructuralDiversity.assign_core_numbers(subgraph)  # O(n+m) Note:n,m of subgraph:
            result[node] = kCoreBaseStructuralDiversity.count_kcores(subgraph, core_numbers, k)  # O(n+m*log(n))
        return result


    class DisjointSet:
        parents: list
        setsize: list

        def __init__(self, length: int):
            """Initialise the parents and setsize arrays."""
            self.parents = [-1] * length
            self.setsize = [1] * length

        def find(self, node: int) -> int:
            """Find the root of the set which contains the given node.
            Use path compression to speed up future searches."""
            if self.parents[node] == -1:
                return node
            else:
                self.parents[node] = self.find(self.parents[node])
                return self.parents[node]

        def union(self, node1: int, node2: int):  # O(log(n))
            """Combine the sets containing node1 and node2. Use union-by-size
            approach, to minimise time complexity."""
            root1 = self.find(node1)
            root2 = self.find(node2)
            if root1 != root2:
                if self.setsize[root1] < self.setsize[root2]:
                    self.setsize[root2] += self.setsize[root1]
                    self.parents[root1] = root2
                else:
                    self.setsize[root1] += self.setsize[root2]
                    self.parents[root2] = root1

    @staticmethod
    def count_kcores(edge_list, core_numbers, k):  # total O(n+m*log(n))
        """count the number of connected components with core number at
        least k using the disjoint set data structure"""
        k_nodes = [n for n in range(len(edge_list)) if core_numbers[n] >= k]  # O(n)
        node_map = {s: i for i, s in enumerate(k_nodes)}  # O(n)

        components = kCoreBaseStructuralDiversity.DisjointSet(len(k_nodes))
        # total for below loops is O(n+m*log(n))
        for node in k_nodes:  # O(n)
            for nbr in edge_list[node]:  # O(deg(node))
                if core_numbers[nbr] >= k:
                    components.union(node_map[node], node_map[nbr])  # O(log(n))

        return len([p for p in components.parents if p == -1])

    @staticmethod
    def get_neighbour_subgraph(edge_list, node):  # total O(n^2)
        """creates a neighbour induced subgraph from node"""
        nbrs = set(edge_list[node])  # O(n)
        node_mapping = {
            nbr: i for i, nbr in enumerate(nbrs)
        }  #O(n)
        subgraph = [set()] * len(nbrs)
        # below loop is bounded by O(n^2)
        for nbr in nbrs:  # O(n)
            subgraph[node_mapping[nbr]] = [
                node_mapping[u] for u in nbrs.intersection(edge_list[nbr])
            ]  # O(min(deg(node), deg(nbr))): bounded theoreticlly by O(n)
        return subgraph

    @staticmethod
    def sort_by_degree(edge_list):  # O(n)
        """initialises the data required for assign_core_numbers, which is the flat array
        implementation of the core decomposition algorithm"""
        n = len(edge_list)
        degrees = [0] * n
        bins = [[] for _ in range(n)]  # nodes can have 0 to n-1 edges
        max_bin = 0
        for node in range(n):  # O(n)
            deg = len(edge_list[node])  # O(1)
            if deg > max_bin:
                max_bin = deg
            degrees[node] = deg
            bins[deg].append(node)

        bin_positions = [0] * (max_bin + 1)
        order = [-1] * n
        node_positions = [-1] * n
        pos = 0
        # both outer and inner loops below are bounded by O(n) for each loop
        # but total number of nodes in bins is n so total is 2n-1
        for deg, nodes in enumerate(bins):  # O(n)
            if deg > max_bin:
                break
            bin_positions[deg] = pos
            for j, node in enumerate(nodes):  # O(n)
                order[pos] = node
                node_positions[node] = pos
                pos += 1

        return order, degrees, bin_positions, node_positions

    @staticmethod
    def assign_core_numbers(edge_list):  # O(n+m)

        """flat array implementation of the core decomposition algorithm"""
        sorted_nodes, node_degrees, degrees_start, node_position =\
             kCoreBaseStructuralDiversity.sort_by_degree(edge_list)  # O(n)

        # total below is O(n+m)
        for i in range(len(sorted_nodes)):  # O(n)
            v = sorted_nodes[i] #the node with the lowest 'degree'
            for u in edge_list[v]: # O(deg(v))
                if node_degrees[u] > node_degrees[v]: #nbr of v with greater degree
                    degree_u = node_degrees[u]
                    position_u = node_position[u]
                    position_w = degrees_start[degree_u]
                    w = sorted_nodes[position_w] #first node with same degree as u
                    if u != w: #swap u and w
                        sorted_nodes[position_u] = w
                        sorted_nodes[position_w] = u
                        node_position[u] = position_w
                        node_position[w] = position_u
                    degrees_start[degree_u] += 1 #shift the boundary where degrees change
                    node_degrees[u] -= 1 #decrement the degree of u
        return node_degrees

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
