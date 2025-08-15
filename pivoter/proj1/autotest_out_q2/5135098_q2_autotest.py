#!/usr/bin/env python3
# Auto-generated for 5135098

STUDENT_ID = "5135098"
STUDENT_NAME = "(Joey) Jingyi Zhang"

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
        
        for v in range(n):
            nbrs = G.adj_list[v]
            if not nbrs:
                sd[v] = 0
                continue

            # 1) 构造邻居诱导子图 H_adj，使用 set 代替 list
            nbr_set = set(nbrs)
            H_adj = {u: set() for u in nbrs}
            for u in nbrs:
                for w in G.adj_list[u]:
                    if w in nbr_set:
                        H_adj[u].add(w)

            # 2) 计算初始度并将 degree < k 的顶点入队
            deg = {u: len(H_adj[u]) for u in nbrs}
            dq = deque([u for u, du in deg.items() if du < k])

            # 3) k-core 剪枝：直接在 nbr_set 中移除度 < k 的顶点
            while dq:
                u = dq.popleft()
                if u not in nbr_set:
                    continue
                nbr_set.remove(u)
                for w in H_adj[u]:
                    if w in nbr_set:
                        deg[w] -= 1
                        if deg[w] < k:
                            dq.append(w)

            # 4) 剩余的 nbr_set 即为 core_vertices
            core_vertices = nbr_set

            # 5) 统计连通分量数
            comp_count = 0
            seen = set()
            for u in core_vertices:
                if u in seen:
                    continue
                comp_count += 1
                stack = [u]
                seen.add(u)
                while stack:
                    x = stack.pop()
                    for w in H_adj[x]:
                        if w in core_vertices and w not in seen:
                            seen.add(w)
                            stack.append(w)

            sd[v] = comp_count

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
