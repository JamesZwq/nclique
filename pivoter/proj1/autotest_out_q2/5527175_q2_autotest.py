#!/usr/bin/env python3
# Auto-generated for 5527175

STUDENT_ID = "5527175"
STUDENT_NAME = "Di Wu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def helper(adj_list):
        num = len(adj_list)
        degs = [len(nb) for nb in adj_list]
        if num == 0:
            return []

        max_deg = max(degs)
        bin = [0] * (max_deg + 1)
        for d in degs:
            bin[d] += 1
        start = 0
        for d in range(max_deg + 1):
            cnt = bin[d]
            bin[d] = start
            start += cnt
        pos = [0] * num
        vert = [0] * num
        for v, d in enumerate(degs):
            idx = bin[d]
            pos[v] = idx
            vert[idx] = v
            bin[d] += 1
        for d in range(max_deg, 0, -1):
            bin[d] = bin[d-1]
        bin[0] = 0
        core = [0] * num
        for i in range(num):
            v = vert[i]
            for u in adj_list[v]:
                if degs[u] > degs[v]:
                    du = degs[u]
                    pu = pos[u]
                    pw = bin[du]
                    w = vert[pw]
                    if u != w:
                        pos[u], pos[w] = pw, pu
                        vert[pu], vert[pw] = w, u
                    bin[du] += 1
                    degs[u] -= 1
            core[v] = degs[v]
        return core

    @staticmethod
    def peel_core(local_deg, k, removal_q, active_flag, epoch, adj):

        while removal_q:
            u = removal_q.popleft()
            if active_flag[u] != epoch:
                continue
            active_flag[u] = 0
            for w in adj[u]:
                if active_flag[w] == epoch:
                    local_deg[w] -= 1
                    if local_deg[w] == k - 1:
                        removal_q.append(w)

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
        n   = G.vertex_num
        adj = G.adj_list
        result  = [0] * n
        if n == 0:
            return []

        if k > 0:
            corenum = kCoreBaseStructuralDiversity.helper(adj)
        else:
            corenum = [0] * n

        active_flag = [0] * n
        visited = [0] * n
        epoch = 1

        for v in range(n):

            nbrs = [u for u in adj[v] if corenum[u] >= k]
            if nbrs:
                for u in nbrs:
                    active_flag[u] = epoch

                local_deg = {
                    u: sum(1 for w in adj[u] if active_flag[w] == epoch)
                    for u in nbrs
                }


                removal_q = deque(u for u, d in local_deg.items() if d < k)

                kCoreBaseStructuralDiversity.peel_core(
                    local_deg, k, removal_q, active_flag, epoch, adj
                )

                complete = 0
                for u in nbrs:
                    if active_flag[u] == epoch and visited[u] != epoch:
                        complete += 1
                        component_q = deque([u])
                        visited[u] = epoch
                        while component_q:
                            x = component_q.popleft()
                            for w in adj[x]:
                                if (active_flag[w] == epoch
                                        and visited[w] != epoch):
                                    visited[w] = epoch
                                    component_q.append(w)
                result[v] = complete

            epoch += 1

        return result



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
