#!/usr/bin/env python3
# Auto-generated for 5574159

STUDENT_ID = "5574159"
STUDENT_NAME = "Xiao Ma"

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
        from collections import deque, defaultdict

        n = G.vertex_num
        adj = G.adj_list

        # Step 1: Compute core numbers for all nodes using peeling
        degree = [len(adj[v]) for v in range(n)]
        core = [0] * n
        removed = [False] * n

        bin_deg = [[] for _ in range(n)]
        for v in range(n):
            bin_deg[degree[v]].append(v)

        curr_deg = 0
        while curr_deg < n:
            while bin_deg[curr_deg]:
                v = bin_deg[curr_deg].pop()
                if removed[v]:
                    continue
                removed[v] = True
                core[v] = curr_deg
                for u in adj[v]:
                    if not removed[u]:
                        d = degree[u]
                        degree[u] -= 1
                        bin_deg[degree[u]].append(u)
            curr_deg += 1

        # Step 2: For each node v, compute τ_k(v)
        result = [0] * n

        for v in range(n):
            neighbors = adj[v]
            nbr_set = set(neighbors)

            # Build neighbor-induced subgraph H = G[N(v)]
            subgraph_adj = defaultdict(list)
            for u in nbr_set:
                for w in adj[u]:
                    if w in nbr_set:
                        subgraph_adj[u].append(w)

            # Step 3: Compute core numbers in neighbor-induced subgraph
            sub_nodes = list(nbr_set)
            sub_idx = {node: i for i, node in enumerate(sub_nodes)}
            sub_n = len(sub_nodes)
            sub_degree = [len(subgraph_adj[node]) for node in sub_nodes]
            sub_core = [0] * sub_n
            sub_removed = [False] * sub_n
            sub_bin_deg = [[] for _ in range(sub_n)]
            for i in range(sub_n):
                sub_bin_deg[sub_degree[i]].append(i)

            curr_d = 0
            while curr_d < sub_n:
                while sub_bin_deg[curr_d]:
                    i = sub_bin_deg[curr_d].pop()
                    if sub_removed[i]:
                        continue
                    sub_removed[i] = True
                    sub_core[i] = curr_d
                    u = sub_nodes[i]
                    for w in subgraph_adj[u]:
                        j = sub_idx[w]
                        if not sub_removed[j]:
                            d = sub_degree[j]
                            sub_degree[j] -= 1
                            sub_bin_deg[sub_degree[j]].append(j)
                curr_d += 1

            # Step 4: Find k-core in neighbor-induced subgraph
            in_k_core = [sub_core[i] >= k for i in range(sub_n)]
            filtered_adj = defaultdict(list)
            for i in range(sub_n):
                if not in_k_core[i]:
                    continue
                u = sub_nodes[i]
                for w in subgraph_adj[u]:
                    j = sub_idx[w]
                    if in_k_core[j]:
                        filtered_adj[u].append(w)

            # Step 5: Count connected components
            visited = set()
            count = 0
            for i in range(sub_n):
                if not in_k_core[i]:
                    continue
                u = sub_nodes[i]
                if u in visited:
                    continue
                count += 1
                queue = deque([u])
                visited.add(u)
                while queue:
                    curr = queue.popleft()
                    for nei in filtered_adj[curr]:
                        if nei not in visited:
                            visited.add(nei)
                            queue.append(nei)

            result[v] = count

        return result

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
