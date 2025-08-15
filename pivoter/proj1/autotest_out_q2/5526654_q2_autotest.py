#!/usr/bin/env python3
# Auto-generated for 5526654

STUDENT_ID = "5526654"
STUDENT_NAME = "Qinhan Xia"

# ======= 学生代码 =======
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        # G.adj_list: List[List[int]]
        # G.vertex_num: int

        n = G.vertex_num
        adj = G.adj_list  # 假设属性名
        tau = [0] * n

        def induced_subgraph(neis):
            # 返回邻居诱导子图的邻接表，节点映射到新编号0..len-1
            neis = list(neis)
            idx_map = {v: i for i, v in enumerate(neis)}
            induced = [[] for _ in neis]
            for i, u in enumerate(neis):
                for v in adj[u]:
                    if v in idx_map and v != u:
                        induced[i].append(idx_map[v])
            return induced

        def find_k_core(adj2, k):
            # 剥壳法求k-core, 返回每个节点是否在k-core里
            n2 = len(adj2)
            deg = [len(neighs) for neighs in adj2]
            alive = [True] * n2
            changed = True
            while changed:
                changed = False
                for i in range(n2):
                    if alive[i] and deg[i] < k:
                        alive[i] = False
                        changed = True
                        for v in adj2[i]:
                            deg[v] -= 1
            return alive

        def count_k_core_components(adj2, k):
            # 返回k-core分量数
            in_core = find_k_core(adj2, k)
            n2 = len(adj2)
            visited = [False] * n2
            count = 0
            for i in range(n2):
                if in_core[i] and not visited[i]:
                    stack = [i]
                    visited[i] = True
                    while stack:
                        u = stack.pop()
                        for v in adj2[u]:
                            if in_core[v] and not visited[v]:
                                visited[v] = True
                                stack.append(v)
                    count += 1
            return count

        for v in range(n):
            neis = set(adj[v])
            if not neis:
                tau[v] = 0
                continue
            adj2 = induced_subgraph(neis)
            tau[v] = count_k_core_components(adj2, k)

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
