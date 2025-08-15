#!/usr/bin/env python3
# Auto-generated for 5511417

STUDENT_ID = "5511417"
STUDENT_NAME = "Dirong Xu"

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
        # TODO
        n = G.vertex_num
        sd = [0] * n

        # 对每个顶点计算其k-core based structural diversity
        for v in range(n):
            # 获取顶点v的邻居
            neighbors = G.adj_list[v] if G.adj_list[v] else []

            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # 构建邻居诱导子图
            neighbor_subgraph = kCoreBaseStructuralDiversity._build_neighbor_subgraph(G, neighbors)

            # 计算k-cores的数量
            sd[v] = kCoreBaseStructuralDiversity._count_k_cores(neighbor_subgraph, k)

        return sd

    # 构建邻居诱导子图
    @staticmethod
    def _build_neighbor_subgraph(G, neighbors):
        """构建由邻居集合诱导的子图"""
        neighbor_set = set(neighbors)
        subgraph = {v: [] for v in neighbors}

        for u in neighbors:
            for w in G.adj_list[u]:
                if w in neighbor_set:
                    subgraph[u].append(w)

        return subgraph

    # 计算k-cores数量
    @staticmethod
    def _count_k_cores(graph, k):
        """计算图中k-cores的数量（连通分量数）"""
        if k == 0:
            # 0-core包含所有顶点，返回连通分量数
            return kCoreBaseStructuralDiversity._count_connected_components(graph)

        # 计算k-core子图
        k_core_graph = kCoreBaseStructuralDiversity._compute_k_core(graph, k)

        if not k_core_graph:
            return 0

        # 计算k-core中的连通分量数
        return kCoreBaseStructuralDiversity._count_connected_components(k_core_graph)

    # 计算k-core子图
    @staticmethod
    def _compute_k_core(graph, k):
        """计算k-core子图，迭代移除度数小于k的顶点"""
        # 复制图结构
        temp_graph = {}
        for v in graph:
            temp_graph[v] = list(graph[v])

        # 迭代移除度数小于k的顶点
        changed = True
        while changed:
            changed = False
            vertices_to_remove = []

            # 找到度数小于k的顶点
            for v in temp_graph:
                if len(temp_graph[v]) < k:
                    vertices_to_remove.append(v)

            # 移除这些顶点
            for v in vertices_to_remove:
                if v in temp_graph:
                    changed = True
                    # 从所有邻居的邻接表中移除v
                    neighbors = temp_graph[v]
                    for u in neighbors:
                        if u in temp_graph and v in temp_graph[u]:
                            temp_graph[u].remove(v)
                    # 移除顶点v
                    del temp_graph[v]

        return temp_graph

    # 计算连通分量数量
    @staticmethod
    def _count_connected_components(graph):
        """使用BFS计算图中连通分量的数量"""
        if not graph:
            return 0

        visited = set()
        component_count = 0

        for start_vertex in graph:
            if start_vertex not in visited:
                # 发现新的连通分量
                component_count += 1

                # BFS遍历整个连通分量
                queue = deque([start_vertex])
                visited.add(start_vertex)

                while queue:
                    current = queue.popleft()
                    for neighbor in graph.get(current, []):
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

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
