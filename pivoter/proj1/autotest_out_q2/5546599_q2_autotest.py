#!/usr/bin/env python3
# Auto-generated for 5546599

STUDENT_ID = "5546599"
STUDENT_NAME = "Yiyu Chen"

# ======= 学生代码 =======
from collections import deque

# v's adj neighbour induced subgraph

def induce_subgraph(G, v, in_k_core):
    # neighbour of v in kcore
    neighbors = [u for u in G.adj_list[v] if in_k_core[u]]
    neighbor_set = set(neighbors) 

    # initializ
    sub_G = {u: [] for u in neighbors}
    for u in neighbors:
        
        sub_G[u] = [nei for nei in G.adj_list[u] if nei in neighbor_set]
    return sub_G

# num of rest conn compo 

def compute_k_core_components(sub_G, k):
    if not sub_G:
        return 0  

    degrees = {u: len(neis) for u, neis in sub_G.items()}
    deleted = set()

    # 
    queue = deque([u for u in sub_G if degrees[u] < k])
    deleted.update(queue)

    
    while queue:
        u = queue.popleft()
        for v in sub_G[u]:
            if v not in deleted:
                degrees[v] -= 1
                if degrees[v] < k:
                    deleted.add(v)
                    queue.append(v)

    # BFS cal the rest conn compo
    visited = set()
    comp_count = 0
    for u in sub_G:
        if u in deleted or u in visited:
            continue
        # found a unvisited node ,start
        comp_count += 1
        q = deque([u])
        visited.add(u)
        while q:
            curr = q.popleft()
            for nei in sub_G[curr]:
                if nei not in deleted and nei not in visited:
                    visited.add(nei)
                    q.append(nei)

    return comp_count


# 
class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
       
        n = G.vertex_num  # num of v
        # degrees
        degree = [len(G.adj_list[v]) for v in range(n)]
        in_k_core = [True] * n  

        # Step 1: top 2 down
        
        changed = True
        while changed:
            changed = False
            for v in range(n):
                if in_k_core[v] and degree[v] < k:
                    # degree 2 small , kcore out
                    in_k_core[v] = False
                    changed = True
                    # delete the adj relationship
                    for u in G.adj_list[v]:
                        if in_k_core[u]:
                            degree[u] -= 1

        # Step 2: 
        tau_k = [0] * n
        for v in range(n):
            if not in_k_core[v]:
                continue  
                
            sub_G = induce_subgraph(G, v, in_k_core)
            tau_k[v] = compute_k_core_components(sub_G, k)

        return tau_k

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
