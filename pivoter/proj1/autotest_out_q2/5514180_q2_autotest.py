#!/usr/bin/env python3
# Auto-generated for 5514180

STUDENT_ID = "5514180"
STUDENT_NAME = "Luze Liu"

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
        ans = [0]*G.vertex_num
        
        for v in range(G.vertex_num): # O(n)
            # 1) collect N(v)
            Nv = set(G.adj_list[v])
            
            # 2) build induced adjacency
            induced = {u: [] for u in Nv}
            for u in Nv:    # O(degm)
                for w in G.adj_list[u]: # O(deg m )
                    if w in Nv:
                        induced[u].append(w)
            
            # 3) peel to k-core
            degree = {u: len(induced[u]) for u in Nv}
            queue = deque(u for u,d in degree.items() if d < k)
            while queue:
                u = queue.popleft()
                for nbr in induced[u]:
                    if degree[nbr] >= k:
                        degree[nbr] -= 1
                        if degree[nbr] == k-1:
                            queue.append(nbr)
                degree[u] = -1  # mark removed
            
            # 4) collect remaining and find connected components
            combine = []
            unvisited = {u for u,d in degree.items() if d >= k}
            queue = deque()
            
            while unvisited:
                node = unvisited.pop()
                box = [node]
                queue.append(node)
                
                while queue:
                    x = queue.popleft()
                    for n in G.adj_list[x]:
                        if n in unvisited:
                            unvisited.remove(n)
                            box.append(n)
                            queue.append(n)
                combine.append(box)
                
            ans[v] = len(combine)
        return ans
                
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
