#!/usr/bin/env python3
# Auto-generated for 5534413

STUDENT_ID = "5534413"
STUDENT_NAME = "Yuankun Jing"

# ======= 学生代码 =======
from collections import deque

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
        sd_list = [0] * n

        # 1. get the core of whole graph
        core_values = kCoreBaseStructuralDiversity.__global_core_decomposition(G.adj_list)

        # 2. Build (k+1) - core subgraph
        important_nodes = []
        for v in range(n):
            if core_values[v] >= k + 1:
                important_nodes.append(v) #core > k+1

        important_set = set()
        for node in important_nodes:
            important_set.add(node)

        # 3. Calculate node by node within the (k+1) - core region
        for v in important_nodes:
            neighbors = []
            for u in G.adj_list[v]:
                if u in important_set:
                    neighbors.append(u)

            subgraph = kCoreBaseStructuralDiversity.__build_subgraph(G.adj_list, neighbors) #Neighbour-induced subgraph
            sd_list[v] = kCoreBaseStructuralDiversity.__compute_local_diversity(subgraph, k)

        return sd_list

    @staticmethod
    def __global_core_decomposition(adj_list):
        """
        Global core number calculation
        """
        n = len(adj_list)  #number of nodes in G

        degree = []
        for v in range(n):
            deg = len(adj_list[v])
            degree.append(deg)

        core = degree.copy()  #core is equal to degree

        bins = []
        for _ in range(n + 1):
          bins.append([]) #Initialize n+1 bins

        for v in range(n):
          bins[degree[v]].append(v)  #Classify points based on degree

        curr_degree = 0
        while curr_degree <= n:
            while bins[curr_degree]:
                u = bins[curr_degree].pop()
                for v in adj_list[u]:
                    if degree[v] > degree[u]:
                        bins[degree[v]].remove(v)  #the neighbor also needs to be affected
                        degree[v] -= 1  #u had been delect
                        bins[degree[v]].append(v)
            curr_degree += 1  #next core bin

        return degree


    @staticmethod
    def __build_subgraph(adj_list, nodes): #adjacency table of the graph
        """
        Construction of neighbor induced subgraph
        """
        node_set = set(nodes)
        subgraph = dict()

        for node in nodes:
            subgraph[node] = set()
            for neighbor in adj_list[node]:
                if neighbor in node_set:
                    subgraph[node].add(neighbor)

        return subgraph

    @staticmethod
    def __compute_local_diversity(subgraph, k):
        """
        Count the number of k-core connected blocks in the current node subgraph
        """
        degrees = dict()
        for node, neighbors in subgraph.items():
            degrees[node] = len(neighbors)
        removed = set()
        queue = deque()
        for node in subgraph:
            if degrees[node] < k:
                queue.append(node)

        while queue:
            node = queue.popleft() #Continuously extract nodes with degrees less than k from the queue and strip them out
            removed.add(node)
            for neighbor in subgraph[node]:
                if neighbor not in removed:
                    degrees[neighbor] -= 1 #After stripping a node, its neighbor degrees need to be reduced by 1
                    if degrees[neighbor] == k - 1:
        #If the neighbor drops to k-1 as a result (i.e. should be stripped off in the next round), add it to the queue
                        queue.append(neighbor)

        active_nodes = [] # Points that satisfy k-core
        for node in subgraph:
            if node not in removed:
                active_nodes.append(node)
        uf = kCoreBaseStructuralDiversity.__UnionFind(active_nodes)

        for node in active_nodes:
            for neighbor in subgraph[node]:
                if neighbor not in removed:
                    uf.union(node, neighbor)

        return uf.count_sets()

    class __UnionFind:
        """
        Local pooling
        """
        def __init__(self, nodes):
          self.parent = {u: u for u in nodes}
          self.rank = {u: 0 for u in nodes}
          self.count = len(nodes)

        def find(self, u): #Find the set where node u is located
          if self.parent[u] != u:
              self.parent[u] = self.find(self.parent[u])
          return self.parent[u]

        def union(self, u, v):#Merge two sets containing u and v
          pu, pv = self.find(u), self.find(v)
          if pu == pv:
            return
          if self.rank[pu] < self.rank[pv]:
            self.parent[pu] = pv
          else:
            self.parent[pv] = pu
            if self.rank[pu] == self.rank[pv]:
              self.rank[pu] += 1
          self.count -= 1

        def count_sets(self): #Return the total number of current collections
            return self.count

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
