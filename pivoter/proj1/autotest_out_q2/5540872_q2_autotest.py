#!/usr/bin/env python3
# Auto-generated for 5540872

STUDENT_ID = "5540872"
STUDENT_NAME = "Yiyao Guo"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute k-core based structural diversity (τ_k) for each vertex using flat-array (CSR) representation.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            Graph with attributes `vertex_num` and `adj_list` (list of neighbor lists).
        k : int
            Core threshold (non-negative).

        Returns
        -------
        List[int]
            τ_k(v) for v in range(n).
        """
        n = G.vertex_num
        if k < 0:
            raise ValueError("k must be non-negative")
        if n == 0:
            return []

        # Build neighbor lists without self-loops or invalids
        nbrs = []
        for u in range(n):
            nbrs.append([v for v in G.adj_list[u] if 0 <= v < n and v != u])

        # Build CSR structure (flat array + offsets)
        offsets = [0] * (n + 1)
        for i in range(n):
            offsets[i + 1] = offsets[i] + len(nbrs[i])
        flat_nbrs = [0] * offsets[n]
        for i in range(n):
            base = offsets[i]
            for j, v in enumerate(nbrs[i]):
                flat_nbrs[base + j] = v

        # Compute core numbers using CSR
        core = kCoreBaseStructuralDiversity._compute_core_numbers_csr(n, offsets, flat_nbrs)

        # Special case k=0: count connected components in neighbor-induced
        if k == 0:
            return [
                kCoreBaseStructuralDiversity._count_components_csr(
                    n, offsets, flat_nbrs, i
                )
                for i in range(n)
            ]

        # Precompute core-neighbors for each v
        core_nbrs = []
        for v in range(n):
            base, end = offsets[v], offsets[v + 1]
            core_nbrs.append([u for u in flat_nbrs[base:end] if core[u] >= k])

        result = [0] * n
        # For each vertex, peel k-core within its neighbor-induced subgraph
        for v in range(n):
            S = set(core_nbrs[v])
            if not S:
                result[v] = 0
                continue

            # Compute degrees within S
            deg = {u: 0 for u in S}
            for u in S:
                bu, eu = offsets[u], offsets[u + 1]
                cnt = 0
                for w in flat_nbrs[bu:eu]:
                    if w in S:
                        cnt += 1
                deg[u] = cnt

            # Peel nodes with degree < k
            q = deque([u for u in S if deg[u] < k])
            removed = set(q)
            while q:
                u = q.popleft()
                bu, eu = offsets[u], offsets[u + 1]
                for w in flat_nbrs[bu:eu]:
                    if w in S and w not in removed:
                        deg[w] -= 1
                        if deg[w] < k:
                            removed.add(w)
                            q.append(w)

            core_nodes = S - removed
            # Count components among core_nodes
            result[v] = kCoreBaseStructuralDiversity._count_components_in_set(
                n, offsets, flat_nbrs, core_nodes
            )
        return result

    @staticmethod
    def _compute_core_numbers_csr(n, offsets, flat_nbrs):
        # Initialize degrees
        deg = [offsets[i + 1] - offsets[i] for i in range(n)]
        if n == 0:
            return []

        # Bin sort
        max_deg = max(deg)
        bin_count = [0] * (max_deg + 1)
        for d in deg:
            bin_count[d] += 1
        start = [0] * (max_deg + 1)
        curr = 0
        for d in range(max_deg + 1):
            start[d] = curr
            curr += bin_count[d]

        # Placement arrays
        pos = [0] * n
        vert = [0] * n
        cnt = [0] * (max_deg + 1)
        for v in range(n):
            d = deg[v]
            idx = start[d] + cnt[d]
            vert[idx] = v
            pos[v] = idx
            cnt[d] += 1

        # Core decomposition
        for idx in range(n):
            v = vert[idx]
            dv = deg[v]
            bv, ev = offsets[v], offsets[v + 1]
            for j in range(bv, ev):
                u = flat_nbrs[j]
                if deg[u] > dv:
                    du = deg[u]
                    pu = pos[u]
                    swap_idx = start[du] + cnt[du] - 1
                    w = vert[swap_idx]
                    if pu < swap_idx:
                        vert[pu], vert[swap_idx] = w, u
                        pos[u], pos[w] = swap_idx, pu
                    cnt[du] -= 1
                    deg[u] -= 1
                    newd = deg[u]
                    cnt[newd] += 1
                    newpos = start[newd] + cnt[newd] - 1
                    vert[newpos] = u
                    pos[u] = newpos
        return deg

    @staticmethod
    def _count_components_csr(n, offsets, flat_nbrs, center):
        # DFS on neighbor-induced subgraph of center
        b, e = offsets[center], offsets[center + 1]
        nbr = flat_nbrs[b:e]
        if not nbr:
            return 0
        visited = set()
        comps = 0
        for u in nbr:
            if u not in visited:
                comps += 1
                stack = [u]
                visited.add(u)
                while stack:
                    w = stack.pop()
                    bw, ew = offsets[w], offsets[w + 1]
                    for x in flat_nbrs[bw:ew]:
                        if x in nbr and x not in visited:
                            visited.add(x)
                            stack.append(x)
        return comps

    @staticmethod
    def _count_components_in_set(n, offsets, flat_nbrs, nodeset):
        if not nodeset:
            return 0
        visited = set()
        comps = 0
        for u in nodeset:
            if u not in visited:
                comps += 1
                stack = [u]
                visited.add(u)
                while stack:
                    w = stack.pop()
                    bw, ew = offsets[w], offsets[w + 1]
                    for x in flat_nbrs[bw:ew]:
                        if x in nodeset and x not in visited:
                            visited.add(x)
                            stack.append(x)
        return comps

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
