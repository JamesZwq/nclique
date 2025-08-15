#!/usr/bin/env python3
# Auto-generated for 5524758

STUDENT_ID = "5524758"
STUDENT_NAME = "Hong Ni"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n
        adj = G.adj_list

        for v in range(n):
            nbrs = adj[v]
            if not nbrs:
                sd[v] = 0
                continue

            # Construct neighbor induction sub-map (including edge deduplication)
            subgraph_nodes = set(nbrs)
            subgraph_adj = defaultdict(list)
            for u in nbrs:
                seen = set() # Prevent repetition of edges
                for w in adj[u]:
                    if w in subgraph_nodes and w not in seen:
                        subgraph_adj[u].append(w)
                        seen.add(w)

            # Calculate k-core
            deg = {u: len(subgraph_adj[u]) for u in subgraph_nodes}
            in_kcore = {u: True for u in subgraph_nodes}
            q = deque([u for u in subgraph_nodes if deg[u] < k])
            while q:
                u = q.popleft()
                if not in_kcore[u]:
                    continue
                in_kcore[u] = False
                for nei in subgraph_adj[u]:
                    if in_kcore[nei]:
                        deg[nei] -= 1
                        if deg[nei] < k:
                            q.append(nei)

            # Find the number of connected components on the stripped sub-picture
            visited = set()
            def bfs(start):
                q = deque([start])
                visited.add(start)
                while q:
                    u = q.popleft()
                    for nei in subgraph_adj[u]:
                        if in_kcore[nei] and nei not in visited:
                            visited.add(nei)
                            q.append(nei)

            count = 0
            for u in subgraph_nodes:
                if in_kcore[u] and u not in visited:
                    bfs(u)
                    count += 1
            sd[v] = count
            
        return sd

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
