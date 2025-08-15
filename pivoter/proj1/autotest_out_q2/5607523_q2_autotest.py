#!/usr/bin/env python3
# Auto-generated for 5607523

STUDENT_ID = "5607523"
STUDENT_NAME = "Yixuan Wang"

# ======= 学生代码 =======
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
        List[int]  # τ_k(v) for all vertices v in G
        """
        n = G.vertex_num
        tau = [0] * n
        # For each vertex v, compute k-core-based structural diversity
        for v in range(n):
            # Build neighbor-induced subgraph on neighbors of v
            nbrs = G.adj_list[v]
            S = set(nbrs)
            # Initialize degree within induced subgraph
            deg = {u: 0 for u in S}
            for u in S:
                # Count neighbors of u within S
                count = 0
                for w in G.adj_list[u]:
                    if w in S:
                        count += 1
                deg[u] = count
            # Remove vertices with degree < k iteratively
            core_nodes = set(S)
            queue = deque([u for u in S if deg[u] < k])
            while queue:
                u = queue.popleft()
                
                if u not in core_nodes:
                    continue
                core_nodes.remove(u)
                # Decrement degree of neighbors
                for w in G.adj_list[u]:
                    if w in core_nodes:
                        deg[w] -= 1
                        if deg[w] < k:
                            queue.append(w)
            # Count connected components in the k-core subgraph
            visited = set()
            count_components = 0
            for u in core_nodes:
                if u not in visited:
                    count_components += 1
                    # BFS from u
                    dq = deque([u])
                    visited.add(u)
                    while dq:
                        x = dq.popleft()
                        for w in G.adj_list[x]:
                            if w in core_nodes and w not in visited:
                                visited.add(w)
                                dq.append(w)
            tau[v] = count_components
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
