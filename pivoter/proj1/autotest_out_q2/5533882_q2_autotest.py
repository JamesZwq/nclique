#!/usr/bin/env python3
# Auto-generated for 5533882

STUDENT_ID = "5533882"
STUDENT_NAME = "Zihan Guo"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from collections import defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def kcore(subgraph, k):

        if len(subgraph) == 0:
            return 0 # No nodes in the subgraph

        deg = defaultdict(dict)
        for u in subgraph:
            deg[u] = len(subgraph[u]) # Get the degree for each node in subgraph
        
        dq = deque()
        dlt = set()

        for v in subgraph:
            if deg[v] < k: # not k-score
                dq.append(v)
                dlt.add(v) # mark as removed

        while dq: # while we have node for the vertex/ it has neightbors that degree less than k
            u = dq.popleft() # pop our from the queue
            for node in subgraph[u]:
                if node in dlt: #check if the node we search is in the dlt list
                    continue
                deg[node] -= 1 # cuz the node is removed, minus the degree of the neighbor
                if deg[node] < k: # after the calculation, if new node' k core is less than k
                    dq.append(node) # push into the queue, continue to delete the node
                    dlt.add(node)

        node_left = []
        for node in subgraph:
            if node not in dlt: # the dlt is finished, the rest node has the degree >= k
                node_left.append(node) # add the node that remain from the degree minus process

        if len(node_left) == 0:
            return 0 # no k-core node left inside the graph

        count = 0
        visited = set()

        for vertex in node_left:
            if vertex not in visited:
                count += 1 # new vertex inside the subgraph, now search the neighbor around it
                q = deque()
                q.append(vertex)
                visited.add(vertex)

                while q:
                    node = q.popleft()# search the neightbor node
                    for nbr in subgraph[node]: # for the left node's neighbor in the subgraph
                        if nbr in node_left and nbr not in visited: # its the new node, not visited
                            visited.add(nbr) # mark as visited, search its next neighbr 
                            q.append(nbr)
        
        return count # as finished the neighbor search, the pq is empty, and get the result of k-core number

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

        for i in range(n):
            nbrs = set(G.adj_list[i])
            if len(nbrs) == 0: # single vertex, no neighbors, deg of 0
                continue # skip this node, look for next node

            subgraph = defaultdict(list)
            for nbr_of_nbr_i in nbrs:
                for node in G.adj_list[nbr_of_nbr_i]: # see whos in the nbr list
                    if node in nbrs: # if the nbr of nbr of i is also the nbr of i
                        subgraph[nbr_of_nbr_i].append(node) # create the subgraph for these node 
            
            sd[i] = kCoreBaseStructuralDiversity.kcore(subgraph, k)

        return sd


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
