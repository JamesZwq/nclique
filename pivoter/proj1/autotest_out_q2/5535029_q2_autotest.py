#!/usr/bin/env python3
# Auto-generated for 5535029

STUDENT_ID = "5535029"
STUDENT_NAME = "Yuan Ou"

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
        n = G.vertex_num
        sd = [0] * n
        #如果图为空，则返回
        for vertex_id,neighbors in enumerate(G.adj_list):
            if not neighbors:
                continue

            #构建仅包含邻居的诱导子图
            neighbor_set = set(neighbors)
            subgraph = {}
            for u in neighbor_set:
                subgraph[u] = []
            for u in neighbor_set:
                for v in G.adj_list[u]:
                    if v in neighbor_set:
                        subgraph[u].append(v)
            #计算 k-core 中连通分量数量
            sd[vertex_id] = kCoreBaseStructuralDiversity.count_core(subgraph, k)
        return sd

    @staticmethod
    def count_core(subgraph, k):
        #如果子图为空，直接返回0
        if not subgraph:
            return 0

        #初始化度数和待删除队列
        degree_set = {}
        for u, neighbor in subgraph.items():
            degree_set[u] = len(neighbor)
        #标记所有度数 < k 的顶点
        remove_core = {u for u, deg in degree_set.items() if deg < k}
        queue = deque(remove_core)

        #只处理被剔除的邻居
        while queue:
            u = queue.popleft()
            for v in (neighbor for neighbor in subgraph[u] if neighbor not in remove_core):
                degree_set[v] -= 1
                if degree_set[v] < k:
                    remove_core.add(v)
                    queue.append(v)

        #剩余顶点列表
        rest_of_core = []
        for u in subgraph:
            if u not in remove_core:
                rest_of_core.append(u)
        #如果没有剩余顶点，返回0
        if not rest_of_core:
            return 0

        #在k-core 中统计连通分量数
        visited_core = set()
        def dfs(cur_core):
            visited_core.add(cur_core)
            for neighbor in subgraph[cur_core]:
                #遍历还没剔除以及未访问的邻居
                if neighbor not in remove_core and neighbor not in visited_core:
                    #往下探索邻居
                    dfs(neighbor)
        num_of_components = 0
        for u in rest_of_core:
            if u not in visited_core:
                #新节点递归调用
                dfs(u)
                num_of_components += 1
        return num_of_components


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
