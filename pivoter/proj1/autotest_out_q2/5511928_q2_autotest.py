#!/usr/bin/env python3
# Auto-generated for 5511928

STUDENT_ID = "5511928"
STUDENT_NAME = "Keyan Miao"

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
        # TODO
        """
          |V|: Number of nodes in the original graph
          |E|: Number of edges in the original graph
          deg(v): Degree of node v
          |V'|: Number of nodes in the neighbor-induced subgraph of v 
          |E'|: Number of edges in the neighbor-induced subgraph of v

        """

        if G.vertex_num == 0:
            return []

        n = G.vertex_num

        # initialize result array
        # time complexity: O(|V|)
        structural_diversity = [0] * n
        
        # total time complexity: O(∑_{v∈V}(deg(v)^2))，i.e., iterate over each node and
        # invoke one O(deg(v)^2) k-core computation per node
        

        # iterate over all nodes and compute the number of k-cores in their neighborhood-induced subgraphs
        # time complexity of the outer loop: O(|V|)
        for v in range(n):

            # get neighbors of node v
            # time complexity: O(1)
            neighbors = G.adj_list[v]
            if not neighbors: 
                structural_diversity[v] = 0
                continue
            
            # the subgraph excludes vertex v itself
            # time complexity: O(deg(v))
            neighbors_set = set(neighbors)
            # time complexity: O(deg(v))
            nbr_induced_sub_G = dict.fromkeys(neighbors_set, [])

            # build neighbor-induced subgraph for node v
            # time complexity: O(deg(v)^2)
            for u in neighbors_set:
              nbr_induced_sub_G[u] = [nbr for nbr in G.adj_list[u] if nbr in neighbors_set]
            
            # call the compute_k_core function to calculate the number of k-cores in the induced subgraph
            # time complexity: O(|V'| + |E'|), approximates to O(deg(v)^2)
            structural_diversity[v] = kCoreBaseStructuralDiversity.compute_k_core(nbr_induced_sub_G, k)
            
        # print(structural_diversity)    
        return structural_diversity

    @staticmethod
    def compute_k_core(nbr_induced_sub_G, k):
        
        if not nbr_induced_sub_G:
            return 0

        # initialize the degree of each node
        # time complexity: O(|E'|), where E' is the number of edges in the neighbor-induced subgraph
        degree = {u: len(adj) for u, adj in nbr_induced_sub_G.items()}

        # initialize the set and queue for nodes to be deleted
        # time complexity: O(|V'|)
        deleted = set()
        queue = deque([u for u in nbr_induced_sub_G if degree[u] < k])
        deleted.update(queue)

        # remove nodes with degree less than k
        # time complexity: O(|E'|)
        while queue:
            u = queue.popleft()
            # O(deg(u))
            for v in nbr_induced_sub_G[u]:
                if v not in deleted:
                    degree[v] -= 1
                    if degree[v] < k:
                        deleted.add(v)
                        queue.append(v)

        # identify remaining nodes that have not been deleted
        # time complexity: O(|V'|)
        valid_nodes = set(nbr_induced_sub_G) - deleted
        if not valid_nodes:
            return 0

        # traverse the remaining nodes to count the number of connected components in the k-core
        # time complexity: O(|V'| + |E'|)
        visited = set()
        count = 0

        for node in valid_nodes:
            if node in visited:
                continue
            count += 1
            queue = deque([node])
            visited.add(node)

            while queue:
                u = queue.popleft()
                # O(deg(u))
                for v in nbr_induced_sub_G[u]:
                    if v not in deleted and v not in visited:
                        visited.add(v)
                        queue.append(v)

        return count
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
