#!/usr/bin/env python3
# Auto-generated for 5497578

STUDENT_ID = "5497578"
STUDENT_NAME = "Shaobo Qiao"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity:
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
            .vertex_num -> number of vertices
            .adj_list   -> adjacency list, list of lists
        k : int
            the k-core threshold

        Returns
        -------
        List[int]: for each v, the number of k-core components in its neighbour subgraph

        Time Complexity
        ---------------
        Let n = |V|, m = |E|.
        For each vertex v:
          - Building the neighbour set: O(deg(v))
          - _peel on H (size h): O(h + sum of edges in H)
          - _count_components on H: O(h + sum of edges in H)
        Summed over all v, total work is O(sum_v deg(v) + sum_v (h_v + e_h_v))
        = O(m + sum_v (h_v + e_h_v)). In sparse graphs, this is roughly O(n + m).
        """
        n = G.vertex_num
        adj = G.adj_list
        result = [0] * n

        def _peel(nodes):
            """Peel off all vertices of degree < k and return the remaining k-core nodes."""
            # Compute initial degrees within the induced subgraph
            deg = {u: sum(1 for w in adj[u] if w in nodes) for u in nodes}
            # Initialize removal queue with vertices whose degree < k
            queue = deque(u for u, d in deg.items() if d < k)
            removed = set()
            while queue:
                u = queue.popleft()
                if u in removed:
                    continue
                removed.add(u)
                # Decrease degree of remaining neighbours and enqueue if they drop below k
                for w in adj[u]:
                    if w in nodes and w not in removed:
                        deg[w] -= 1
                        if deg[w] == k - 1:
                            queue.append(w)
            # Remaining nodes form the k-core
            return nodes - removed

        def _count_components(nodes):
            """Count how many connected components are in the given set of nodes."""
            seen = set()
            comps = 0
            for u in nodes:
                if u not in seen:
                    comps += 1
                    dq = deque([u])
                    seen.add(u)
                    # BFS to mark all nodes in this component
                    while dq:
                        x = dq.popleft()
                        for w in adj[x]:
                            if w in nodes and w not in seen:
                                seen.add(w)
                                dq.append(w)
            return comps

        # Main loop: for each vertex, process its neighbour-induced subgraph
        for v in range(n):
            nbrs = set(adj[v])  # neighbour set H
            if not nbrs:
                # No neighbours => diversity is 0
                continue
            core_nodes = _peel(nbrs)
            # Count how many k-core components remain
            result[v] = _count_components(core_nodes)

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
