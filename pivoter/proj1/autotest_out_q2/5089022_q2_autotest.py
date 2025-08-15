#!/usr/bin/env python3
# Auto-generated for 5089022

STUDENT_ID = "5089022"
STUDENT_NAME = "Quan Zhang"

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
        # TODO

        # Fucntion to build subgraph for selected node and its nbrs 
        def build_graph(nodes, G):
            map_idx = {u:i for i,u in enumerate(nodes)}
            degree = len(nodes)

            sub_degree = [0]*degree # initialise degrees of subgraph as zeros
            adj_list = [[] for _ in range(degree)] # inisitlise adj list for subgraph

            # build subgraph adjacency and degrees
            for i,u in enumerate(nodes):
                # Scan nodes in nbrs_set, if its nbr is also in subgraph, then update the list
                for w in G.adj_list[u]:
                    j = map_idx.get(w, -1) # if not insde the subgraph, j = -1
                    if j >= 0:
                        adj_list[i].append(j)
                        sub_degree[i] += 1
        
            return adj_list, sub_degree
        
        n = G.vertex_num   
        sd = [0]*n   
        
        for v in range(n):
            # get neighbours of v, store as set (unique)
            nbrs_set = set(G.adj_list[v])

            if not nbrs_set:   # no nbr nodes
                sd[v] = 0
                continue

            degree = len(nbrs_set) # number of nbrs
            if k > 0 and degree < k:     # cannot form any k-core
                sd[v] = 0
                continue

            nbrs = list(nbrs_set)
            adj, sub_degree = build_graph(nbrs,G)

            # initialise removed list, when marked as True means not the node forming k-core
            removed = [0]*degree    

            # if a node inside subgraph has degree < k, add it to queue in initialisation
            q = deque(i for i in range(degree) if sub_degree[i] < k)
            while q:
                u = q.popleft()
                if removed[u]:
                    continue
                removed[u] = 1
                for w in adj[u]:   # update nbrs
                    if not removed[w]:
                        sub_degree[w] -= 1    # update nbrs degree
                        if sub_degree[w] < k:    # if after update, degree falls under k, add it in the quete to be removed
                            q.append(w)

            # update on sd[v] for the final counting
            count = 0
            # for subgraph bfs
            visited = [0]*degree

            for i in range(degree):
                # removed nodes skipped, so not incrementing count
                if not removed[i] and not visited[i]:
                    count += 1
                    stack = [i]
                    visited[i] = 1
                    while stack:
                        x = stack.pop()
                        for y in adj[x]:
                            if not removed[y] and not visited[y]:
                                visited[y] = 1
                                stack.append(y)

            sd[v] = count

        return sd

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
