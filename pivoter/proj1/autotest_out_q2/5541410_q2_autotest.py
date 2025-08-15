#!/usr/bin/env python3
# Auto-generated for 5541410

STUDENT_ID = "5541410"
STUDENT_NAME = "Linjun Pang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
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

        # generate a subgraph of each vertex's neighbours
        for v in range(n):
            nbrs = set(G.adj_list[v])
            # if no neighbors, continue for better time complexity
            if len(nbrs) == 0:
                continue
            
            # create hashmap+hashset to store the adjacency list of the subgraph
            subgraph = defaultdict(set)
            # add the edges between v's the neighbours
            for u in nbrs:
                for w in G.adj_list[u]:
                    if w in nbrs and w != u:
                        subgraph[u].add(w)

            # count the tau of k-core for each vertex
            sd[v] = kCoreBaseStructuralDiversity.count_kcore_tau(subgraph, k)

        return sd
    
    def count_kcore_tau(subgraph, k):
        vertexs = set(subgraph.keys())
        removed_vertexs = set()

        # count the degree of each vertex in subgraph
        subgraph_degree = {}
        for v in vertexs:
            subgraph_degree[v] = len(subgraph[v])

        # put the vertexs whose degree is less than k in a list for further process
        processing_nodes = []
        for v in vertexs:
            if subgraph_degree[v] < k:
                processing_nodes.append(v)

        while processing_nodes:
            # create a list to store the nodes which will be removed in the next loop
            to_remove = []

            # processing the nodes that degree less than k
            for v in processing_nodes:
                if v in removed_vertexs:
                    continue

                removed_vertexs.add(v)

                for nbr in subgraph[v]:
                    # each of it's nbr's degree will be -1
                    if nbr not in removed_vertexs:
                        subgraph_degree[nbr] -= 1
                        # if a nbr node's degree is less than k after current node is removed, put it in the list
                        if subgraph_degree[nbr] < k:
                            to_remove.append(nbr)
            # update the process list
            processing_nodes = to_remove

        # store the rest vertexs for counting the cc
        final_vertexs = vertexs - removed_vertexs
        if not final_vertexs:
            return 0
        
        return kCoreBaseStructuralDiversity.count_connected_components(subgraph, final_vertexs)
    
    def count_connected_components(subgraph, final_vertexs):
        # create indexs for each vertexs for union operations
        vertex_ids = {}
        for i, v in enumerate(final_vertexs):
            vertex_ids[v] = i
        n = len(final_vertexs)

        # union-find algorithm to find all the connected components
        parent = [i for i in range(n)]
        size   = [1 for i in range(n)]
        def find(e):
            if parent[e] == e:
                return e
            else:
                parent[e] = find(parent[e]) 
                return parent[e]
        
        def union(x, y):
            x = find(x)
            y = find(y)
            if x == y:
                return
            if size[x] < size[y]:
                size[y] += size[x]
                parent[x] = y
            else:
                size[x] += size[y]
                parent[y] = x

        # create connected component
        for u in final_vertexs:
            u_idx = vertex_ids[u]
            for v in subgraph[u]:
                # avoid duplicate union
                if v in final_vertexs and u < v:
                    v_idx = vertex_ids[v]
                    union(u_idx, v_idx)
        
        # create a set to store the roots
        roots = set()
        for i in range(n):
            roots.add(find(i))
        
        # the number of root is the number of cc
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
