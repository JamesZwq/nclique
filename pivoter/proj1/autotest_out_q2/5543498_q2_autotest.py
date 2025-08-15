#!/usr/bin/env python3
# Auto-generated for 5543498

STUDENT_ID = "5543498"
STUDENT_NAME = "Jie Sha"

# ======= 学生代码 =======
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        total_nodes = G.vertex_num
        diversity_scores = [0] * total_nodes
        for node in range(total_nodes):
            neighbors = G.adj_list[node]
            if not neighbors:
                continue
            subgraph_nodes = list(set(neighbors))
            if not subgraph_nodes:
                continue
            subgraph_adj, subgraph_deg = kCoreBaseStructuralDiversity._construct_subgraph(G.adj_list, subgraph_nodes)
            valid_nodes = kCoreBaseStructuralDiversity._apply_kcore_pruning(subgraph_adj, subgraph_deg, k)
            connected_parts = kCoreBaseStructuralDiversity._explore_components(subgraph_adj, valid_nodes)
            diversity_scores[node] = connected_parts
        return diversity_scores

    @staticmethod
    def _construct_subgraph(global_adj, vertex_list):
        idx_map = {v: i for i, v in enumerate(vertex_list)}
        size = len(vertex_list)
        sub_adj = [[] for _ in range(size)]
        sub_deg = [0] * size
        for i, u in enumerate(vertex_list):
            for v in global_adj[u]:
                if v in idx_map:
                    j = idx_map[v]
                    sub_adj[i].append(j)
                    sub_deg[i] += 1
        return sub_adj, sub_deg

    @staticmethod
    def _apply_kcore_pruning(adj, degrees, k):
        n = len(adj)
        alive = [True] * n
        pending = deque(i for i in range(n) if degrees[i] < k)
        for i in pending:
            alive[i] = False
        while pending:
            current = pending.popleft()
            for neighbor in adj[current]:
                if alive[neighbor]:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        alive[neighbor] = False
                        pending.append(neighbor)
        return alive

    @staticmethod
    def _explore_components(adj, active):
        n = len(adj)
        seen = [False] * n
        component_count = 0
        for i in range(n):
            if active[i] and not seen[i]:
                component_count += 1
                queue = deque([i])
                seen[i] = True
                while queue:
                    u = queue.popleft()
                    for v in adj[u]:
                        if active[v] and not seen[v]:
                            seen[v] = True
                            queue.append(v)
        return component_count
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
