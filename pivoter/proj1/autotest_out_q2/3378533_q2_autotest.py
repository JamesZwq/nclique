#!/usr/bin/env python3
# Auto-generated for 3378533

STUDENT_ID = "3378533"
STUDENT_NAME = "Jian Song"

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
        n = G.vertex_num
        sd = [0] * n

        # we loop through each vertex v to calculate its tau value
        for v in range(n):

            # we store vertex v's neighbors in a hash table for O(1) lookup
            v_adj = set(G.adj_list[v])
            
            # we first check if vertex v has any neighbors, if not, then its tau remains 0
            if not v_adj:
                continue

            # we build sub-graph with only neighbors of vertex v
            H_adj = {}
            for u in v_adj:
                H_adj[u] = []
                
                for w in G.adj_list[u]:
                    if w in v_adj:
                        H_adj[u].append(w)


            # we then calculate k-core of the sub-graph using a queue structure
            # First, we calculate degree of each vertex in the sub-graph and then add the vertices whose degree is less than k to the queue
            H_deg = {}
            queue = deque()
            for u in v_adj:
                deg_u = len(H_adj[u])
                H_deg[u] = deg_u
                
                if deg_u < k:
                    queue.append(u)

            # We will now remove the vertices whose degree is less than k
            # we store removed vertices in a hash table for O(1) lookup
            removed = set()
            while queue:
                u = queue.popleft()
                if u in removed:
                    continue
                removed.add(u)

                # update degree of each vertex in sub-graph H after vertex u is removed
                for w in H_adj[u]:
                    if w not in removed:
                        H_deg[w] -= 1
                        if H_deg[w] == k - 1:
                            queue.append(w)

            # we now update sub-graph of vertex v without vertices whose degree is less than k
            v_kcore = []
            for u in v_adj:
                if u not in removed:
                    v_kcore.append(u)
            
            # if the sub-graph is empty, then vertex v's k-core value is 0
            if not v_kcore:
                continue

            # we now run DFS using a stack (python list) to count the connected components
            visited = set()
            conn_comp = 0
            for u in v_kcore:
                if u not in visited:
                    conn_comp += 1
                    
                    stack = [u]
                    visited.add(u)
                    while stack:
                        x = stack.pop()
                        for w in H_adj[x]:
                            if w not in removed and w not in visited:
                                visited.add(w)
                                stack.append(w)
            sd[v] = conn_comp

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
