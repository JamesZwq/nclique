#!/usr/bin/env python3
# Auto-generated for 5528464

STUDENT_ID = "5528464"
STUDENT_NAME = "Jingjing Xia"

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
        tau = [0] * G.vertex_num

        for v in range(G.vertex_num):
            neighbors = G.adj_list[v]

            if not neighbors:
                tau[v] = 0
                continue
            subgraph = defaultdict(list)
            node_set = set(neighbors)
            node_mapping = {u: i for i, u in enumerate(neighbors)}

            for u in neighbors:
                for neighbor in G.adj_list[u]:
                    if neighbor in node_set and neighbor != v:
                        subgraph[node_mapping[u]].append(node_mapping[neighbor])

            if not subgraph:
                tau[v] = 0
                continue

            n = len(neighbors)
            degrees = [0] * n
            for u in subgraph:
                degrees[u] = len(subgraph[u])

            max_degree = max(degrees) if degrees else 0
            bins = [0] * (max_degree + 1)
            for d in degrees:
                bins[d] += 1

            start = 0
            for d in range(max_degree + 1):
                num = bins[d]
                bins[d] = start
                start += num

            pos = [0] * n
            vert = [0] * n
            for u in range(n):
                pos[u] = bins[degrees[u]]
                vert[pos[u]] = u
                bins[degrees[u]] += 1

            for d in range(max_degree, 0, -1):
                bins[d] = bins[d - 1]
            bins[0] = 0

            core = [0] * n
            for i in range(n):
                u = vert[i]
                core[u] = degrees[u]
                for v_neighbor in subgraph.get(u, []):
                    if degrees[v_neighbor] > degrees[u]:
                        dv = degrees[v_neighbor]
                        pos_v = bins[dv]
                        neighbor_pos = pos[v_neighbor]

                        if pos_v != neighbor_pos:
                            vert[neighbor_pos], vert[pos_v] = vert[pos_v], vert[neighbor_pos]
                            pos[vert[neighbor_pos]] = neighbor_pos
                            pos[vert[pos_v]] = pos_v

                        bins[dv] += 1
                        degrees[v_neighbor] -= 1

            visited = [False] * n
            core_count = 0

            for u in range(n):
                if not visited[u] and core[u] >= k:
                    queue = deque()
                    queue.append(u)
                    visited[u] = True
                    core_count += 1

                    while queue:
                        current = queue.popleft()
                        for neighbor in subgraph.get(current, []):
                            if not visited[neighbor] and core[neighbor] >= k:
                                visited[neighbor] = True
                                queue.append(neighbor)

            tau[v] = core_count

        return tau

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
