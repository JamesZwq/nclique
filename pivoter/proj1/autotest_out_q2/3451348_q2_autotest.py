#!/usr/bin/env python3
# Auto-generated for 3451348

STUDENT_ID = "3451348"
STUDENT_NAME = "Wencong Zuo"

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

        # If graph is empty, return empty tau
        if G.vertex_num == 0:
            return []

        def extractNeighbour(G, v):
            """
                Extract sub-graph of v, which only contains the neighbour nodes and current connection among them.
                Args:
                    G (UndirectedUnweightedGraph instance): The graph need to be searched
                    v (int): id of specific node
            """
            # get the neighbour of v
            Neighbour = set(G.adj_list[v])
            # flags of neighbour list
            is_neighbour = [False] * G.vertex_num
            for u in Neighbour:
                is_neighbour[u] = True

            # generate sub-graph edge list with meta_data as first element.
            H_edge_list = [[len(Neighbour), 0]]
            for u in Neighbour:
                for w in G.adj_list[u]:
                    if is_neighbour[w] and u < w:
                        H_edge_list.append([u, w])
            # update metatdata
            H_edge_list[0][1] = len(H_edge_list) - 1
            return H_edge_list
            # Time Complexity O(sum(d(u)) for u in Neighbour)

        def dfs(u, adjacency, current):
            """
                function to find connected nodes for u, with restriction of currrent node list
                Args:
                    u (int): id for the original node
                    adjacency (dict): adjacency list
                    current (set): current node list
            """
            # initialise stack and visited set
            stack = [u]
            visited = {u}
            # iterate until empty stack
            while stack:
                # pop element from stack as current
                x = stack.pop()
                # find out next node that not visited
                for w in adjacency[x]:
                    if w in current and w not in visited:
                        visited.add(w)
                        stack.append(w)
            return visited
            # Time complexity O(d(v) + m(v))
        sd = [] # initialise tau list

        # iterate through all nodes
        for v in range(G.vertex_num):
            # tau = 0 if no adjacency nodes for v
            if G.adj_list[v] == []:
                sd.append(0)
                continue

            # find out the sub_graph H
            H = extractNeighbour(G, v)

            # generate adjacency list for H
            H_adj_list = {}
            for u, v in H[1:]:
                H_adj_list[u] = H_adj_list.get(u, []) + [v]
                H_adj_list[v] = H_adj_list.get(v, []) + [u]

            # generate degree list for H
            H_degree_list = {}
            for u in H_adj_list:
                H_degree_list[u] = len(H_adj_list[u])

            # if k is larger than max degree for H, then tau = 0
            if k > max(H_degree_list.values()):
                sd.append(0)
                continue

            # initialise a queue and get all nodes in H
            queue = deque()
            current_node = set(H_adj_list.keys())
            # drop all nodes with degree smaller than k to queue
            for u in H_degree_list:
                if H_degree_list[u] < k:
                    queue.append(u)
            # iterate until empty queue
            while queue:
                # remove u from current node set if not yet and reduce the degree for
                # node in current node
                u = queue.pop()
                if u not in current_node:
                    continue
                current_node.remove(u)
                for v in H_adj_list[u]:
                    if v in current_node:
                        H_degree_list[v] -= 1
                        if H_degree_list[v] < k:
                            queue.append(v)
            # run DFS, with limitation of current node set
            visited = set() # initialise visited set
            component_count = 0 # initialise count
            # iterate a node on each connected components, and count how many dfs have
            # been done.
            for u in current_node:
                if u not in visited:
                    component_count += 1
                    visited = visited | dfs(u, H_adj_list, current_node)
            sd.append(component_count)
        return sd
        # Time complexity O(sum(d(u)^2) where u belongs to vertex)

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
