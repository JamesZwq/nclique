#!/usr/bin/env python3
# Auto-generated for 5507302

STUDENT_ID = "5507302"
STUDENT_NAME = "Yiyang Xu"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        tau = [0]*n

        def count_k_cores_in_subgraph(adj, k):
            # input: adj - induced adjacency list (neigh-indexed), k
            node_num = len(adj)
            if node_num == 0:
                return 0

            remaining = set(range(node_num))
            count = 0

            while remaining:
                # Build current subgraph node id mapping (for this round)
                old2new = {x: i for i, x in enumerate(remaining)}
                new2old = list(remaining)
                subadj = [[] for _ in range(len(remaining))]
                for ni, oi in enumerate(new2old):
                    for oj in adj[oi]:
                        if oj in remaining:
                            subadj[ni].append(old2new[oj])

                # K-core peeling in this subgraph
                deg = [len(neigh) for neigh in subadj]
                alive = [True]*len(deg)
                q = deque([i for i, d in enumerate(deg) if d < k])
                while q:
                    u = q.popleft()
                    alive[u] = False
                    for v in subadj[u]:
                        if alive[v]:
                            deg[v] -= 1
                            if deg[v] == k-1:
                                q.append(v)
                # Find all connected components among alive nodes
                vis = [False]*len(alive)
                for i in range(len(alive)):
                    if alive[i] and not vis[i]:
                        # BFS to get the component
                        comp = []
                        dq = deque([i])
                        vis[i] = True
                        while dq:
                            u = dq.popleft()
                            comp.append(u)
                            for v in subadj[u]:
                                if alive[v] and not vis[v]:
                                    vis[v] = True
                                    dq.append(v)
                        # For this component, check if all nodes have deg ≥ k in the *component-induced* subgraph
                        # Build comp-subadj and peel in this component, recursively if necessary!
                        if len(comp) > 0:
                            comp_size = len(comp)
                            cidx = {x: idx for idx,x in enumerate(comp)}
                            compadj = [[] for _ in range(comp_size)]
                            for ii, ni in enumerate(comp):
                                for nj in subadj[ni]:
                                    if nj in comp:
                                        compadj[ii].append(cidx[nj])
                            # check local degs
                            cdeg = [len(lst) for lst in compadj]
                            if min(cdeg) >= k:
                                count += 1
                                for ni in comp:
                                    remaining.discard(new2old[ni])
                            else:
                                # Need to do recursive peeling; so replace remaining with this component & break to restart on it
                                remaining = set([new2old[ni] for ni in comp])
                                break
                else:
                    # All components processed this round; break outer loop
                    break
            return count

        for v in range(n):
            nbrs = G.adj_list[v]
            m = len(nbrs)
            if m == 0:
                tau[v] = 0
                continue
            idmap = {u:i for i,u in enumerate(nbrs)}
            nadj = [[] for _ in range(m)]
            for i,u in enumerate(nbrs):
                for w in G.adj_list[u]:
                    j = idmap.get(w)
                    if j is not None:
                        nadj[i].append(j)
            tau[v] = count_k_cores_in_subgraph(nadj, k)
        return tau

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
