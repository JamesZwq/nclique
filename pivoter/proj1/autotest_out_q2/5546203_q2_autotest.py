#!/usr/bin/env python3
# Auto-generated for 5546203

STUDENT_ID = "5546203"
STUDENT_NAME = "Jinghuang Zhang"

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
        n = G.vertex_num
        sd = [0] * n
        core = kCoreBaseStructuralDiversity._coreNumber(G)
        
        mark = [False] * n
        alive = [False] * n
        deg = [0] * n
        
        for v in range(n): 
            neighbors = G.adj_list[v]
            if len(neighbors) < k:
                continue
            
            candidates = []
            for u in neighbors:
                if core[u] >= k:    # Skip if the core of u is smaller than k in the original graph
                    if not mark[u]:
                        mark[u] = True
                        candidates.append(u)
                        
            if len(candidates) < k: # Skip if the number of nodes in candidates is smaller than k
                for u in candidates:
                    mark[u] = False
                continue
            
            peel = deque()
            for u in candidates:
                c = 0
                for w in G.adj_list[u]:
                    if mark[w]:
                        c += 1
                deg[u] = c
                alive[u] = True
                if c < k:
                    peel.append(u)
            
            while peel:
                u = peel.popleft()
                if not alive[u]:
                    continue
                alive[u] = False
                for w in G.adj_list[u]:
                    if alive[w]:
                        deg[w] -= 1
                        if deg[w] == k - 1:
                            peel.append(w)
            
            comps = 0
            bfs = deque()
            for u in candidates:
                if alive[u]:
                    comps += 1
                    alive[u] = False
                    bfs.append(u)
                    while bfs:
                        x = bfs.popleft()
                        for w in G.adj_list[x]:
                            if alive[w]:
                                alive[w] = False
                                bfs.append(w)
            
            sd[v] = comps
        
            for u in candidates:
                mark[u] = False
                deg[u] = 0
                
        return sd

    @staticmethod
    def _coreNumber(G):
        n = G.vertex_num
        degree = [len(G.adj_list[u]) for u in range(n)]
        maxDegree = max(degree) if n else 0
        bins = [0 for _ in range(maxDegree + 1)]
        for d in degree:
            bins[d] += 1
        start = 0
        for d in range(maxDegree+1):
            c = bins[d]
            bins[d] = start
            start += c 
        pos = [0 for _ in range(n)] # the position of node u in vert array
        vert = [0 for _ in range(n)] 
        for u in range(n):
            d = degree[u]
            pos[u] = bins[d]
            vert[pos[u]] = u
            bins[d] += 1
        for d in range(maxDegree,0,-1):
            bins[d] = bins[d-1]
        bins[0] = 0

        core = [0]*n
        for i in range(n):
            u = vert[i]
            uDegree = degree[u]
            core[u] = uDegree
            for w in G.adj_list[u]:
                if degree[w] > uDegree:
                    wDegree = degree[w]
                    pDegree = pos[w]
                    pDegree_new = bins[wDegree]
                    if pDegree != pDegree_new:
                        vw = vert[pDegree_new]
                        vert[pDegree_new], vert[pDegree] = vert[pDegree], vert[pDegree_new]
                        pos[w] = pDegree_new
                        pos[vw] = pDegree
                    bins[wDegree] += 1
                    degree[w] -= 1
        return core

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
