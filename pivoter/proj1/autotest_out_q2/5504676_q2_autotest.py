#!/usr/bin/env python3
# Auto-generated for 5504676

STUDENT_ID = "5504676"
STUDENT_NAME = "Yao Zhu"

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
        adj = G.adj_list
        # 1. globally calculate the Core Number of each vertex
        core = kCoreBaseStructuralDiversity._compute_core_numbers(adj)
        # initialization results
        tau = [0] * n
        # 2. find all vertices that satisfy the global core number >= k+1, and construct the vertex set of the (k+1)-core subgraph
        threshold = k + 1
        core_nodes = {v for v, c in enumerate(core) if c >= threshold}
        if not core_nodes:
            return tau
        # 3. find all connected components in the (k+1)-core subgraph
        components = kCoreBaseStructuralDiversity._find_components(core_nodes, adj)
        # 4. for each connected component batch, update tau
        for comp in components:
            kCoreBaseStructuralDiversity._process_component_batch(comp, k, core_nodes, adj, tau)
        return tau

    def _compute_core_numbers(adj):
        ### compute core number for each vertex via bucket algorithm in O(n+m)
        n = len(adj)
        deg = [len(adj[v]) for v in range(n)]
        if n == 0:
            return []
        maxd = max(deg)
        # count sort by degree
        bin_count = [0] * (maxd + 1)
        for d in deg:
            bin_count[d] += 1
        start = [0] * (maxd + 1)
        s = 0
        for d in range(maxd + 1):
            start[d] = s
            s += bin_count[d]
        vert = [0] * n
        pos = [0] * n
        # place vertices
        nxt = start.copy()
        for v in range(n):
            d = deg[v]
            pos[v] = nxt[d]
            vert[nxt[d]] = v
            nxt[d] += 1
        core = [0] * n
        # peel
        for i in range(n):
            v = vert[i]
            core[v] = deg[v]
            for u in adj[v]:
                if deg[u] > deg[v]:
                    du = deg[u]
                    pu = pos[u]
                    pw = start[du]
                    w = vert[pw]
                    if u != w:
                        # swap u and w in vert
                        pos[u], pos[w] = pw, pu
                        vert[pu], vert[pw] = w, u
                    start[du] += 1
                    deg[u] -= 1
        return core

    def _find_components(nodes, adj):
        ### find connected components within 'nodes' subset
        visited = set()
        comps = []
        for v in nodes:
            if v in visited:
                continue
            dq = deque([v])
            visited.add(v)
            comp = []
            while dq:
                u = dq.popleft()
                comp.append(u)
                for w in adj[u]:
                    if w in nodes and w not in visited:
                        visited.add(w)
                        dq.append(w)
            comps.append(comp)
        return comps

    def _process_component_batch(comp, k, core_nodes, adj, tau):
        ### for a given (k+1)-core connected component, update tau[v] = max(tau[v], local structural diversity from this comp). Process both vertices in comp and their neighbors.
        comp_set = set(comp)
        # related vertices: those within comp and those adjacent to comp
        related = set(comp)
        for u in comp:
            for v in adj[u]:
                related.add(v)
        # calculate structural diversity for each relevant vertex
        for v in related:
            # only consider the intersection of v's neighbors with the current comp
            nbrs = [w for w in adj[v] if w in comp_set]
            if not nbrs:
                sd = 0
            else:
                # local k-core peeling
                sd = kCoreBaseStructuralDiversity._compute_k_core_count(nbrs, k, adj)
            # take the maximum value
            if sd > tau[v]:
                tau[v] = sd

    def _compute_k_core_count(nodes, k, adj):
        ### given a list of nodes (subgraph), peel k-core and count its connected pieces
        S = set(nodes)
        if len(S) == 0:
            return 0
        # constructing degree and sub-adjacency lists
        deg = {}
        sub_adj = {}
        for u in S:
            nbrs = [w for w in adj[u] if w in S]
            sub_adj[u] = nbrs
            deg[u] = len(nbrs)
        # strip
        queue = deque([u for u in S if deg[u] < k])
        removed = set(queue)
        while queue:
            u = queue.popleft()
            for w in sub_adj[u]:
                if w not in removed:
                    deg[w] -= 1
                    if deg[w] < k:
                        removed.add(w)
                        queue.append(w)
        core_nodes = S - removed
        # counting connected components
        visited = set()
        count = 0
        for u in core_nodes:
            if u not in visited:
                count += 1
                dq = deque([u])
                visited.add(u)
                while dq:
                    x = dq.popleft()
                    for w in sub_adj[x]:
                        if w in core_nodes and w not in visited:
                            visited.add(w)
                            dq.append(w)
        return count


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
