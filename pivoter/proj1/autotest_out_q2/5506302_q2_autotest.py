#!/usr/bin/env python3
# Auto-generated for 5506302

STUDENT_ID = "5506302"
STUDENT_NAME = "Wenqiang Tan"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################


class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def neighbor_induced_subgraph(adj_list, vertex):
        nbrs = adj_list[vertex]
        subgraph = {}
        for i in range(len(adj_list)):
            if i == vertex or i not in nbrs:
                continue
            temp = []
            for v in adj_list[i]:
                if v in nbrs:
                    temp.append(v)
            if i in nbrs:
                subgraph[i] = temp
        return subgraph

    @staticmethod
    def k_core_decomposition(adj_list, k):
        # recursively remove vertices that degree less than k
        changed = True
        while changed:
            changed = False
            remove = [node for node, nbrs in adj_list.items() if len(nbrs) < k]
            if remove:
                changed = True
                for node in remove:
                    del adj_list[node]
                for node in adj_list:
                    adj_list[node] = [nbr for nbr in adj_list[node] if nbr not in remove]
        return adj_list
    
    @staticmethod
    def connected_components(subgraph):
        visited = set()
        components = []
        def dfs(node, comp):
            visited.add(node)
            comp.append(node)
            for nbr in subgraph[node]:
                if nbr not in visited:
                    dfs(nbr, comp)
        for node in subgraph:
            if node not in visited:
                comp = []
                dfs(node, comp)
                components.append(comp)
        return components

    @staticmethod
    def tau(adj_list, v, k):
        subgraph = kCoreBaseStructuralDiversity.neighbor_induced_subgraph(adj_list, v)
        subgraph = kCoreBaseStructuralDiversity.k_core_decomposition(dict(subgraph), k)
        cc = kCoreBaseStructuralDiversity.connected_components(subgraph)
        return len(cc)

    @staticmethod
    def process(G, k):
        n = G.vertex_num
        sd = [kCoreBaseStructuralDiversity.tau(G.adj_list, i, k) for i in range(n)]
        return sd


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
