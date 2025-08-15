#!/usr/bin/env python3
# Auto-generated for 5514965

STUDENT_ID = "5514965"
STUDENT_NAME = "Dinan Jiang"

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

        n = G.vertex_num  # 图中节点数 Number of vertices
        sd = [0] * n      # 初始化返回数组，每个位置是一个节点的结构多样性 Initialize result list

        for v in range(n):  # 遍历每个节点 Loop over each node
            neighbors = G.adj_list[v]  # 获取当前节点的邻居 Get neighbors of node
            if not neighbors:
                continue  # 如果没有邻居，跳过 No neighbors, τ_k(v)=0

            # 构造邻居诱导子图 Build neighbor-induced subgraph
            subgraph = {u: [] for u in neighbors}
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in subgraph:
                        subgraph[u].append(w)

            # 对子图进行k-core剥离 Extract k-core from subgraph
            core = kCoreBaseStructuralDiversity.compute_k_core(subgraph, k)

            # 统计连通块数量 Count number of connected components
            sd[v] = kCoreBaseStructuralDiversity.count_connected_components(core)

        return sd  # 返回结构多样性列表 Return the result list

    @staticmethod
    def compute_k_core(graph, k):
        """
        输入子图，输出它的k-core
        Remove nodes with degree < k until stable
        """
        degrees = {u: len(graph[u]) for u in graph}  # 初始度数 Initialize degrees
        queue = deque([u for u in graph if degrees[u] < k])  # 队列：要剥离的点 Queue of nodes to remove

        while queue:
            u = queue.popleft()  # 取出要删除的点 Pop node to remove
            for v in graph[u]:  # 遍历它的邻居 Visit its neighbors
                if v in graph:
                    if u in graph[v]:
                        graph[v].remove(u)  # 更新邻居的邻接表 Remove u from neighbor
                        degrees[v] -= 1     # 减少度数 Decrease degree
                        if degrees[v] == k - 1:
                            queue.append(v)  # 新的低度点加入队列 Add to queue
            del graph[u]  # 删除当前点 Remove u from graph

        return graph  # 返回剩下的子图 Return the final k-core subgraph

    @staticmethod
    def count_connected_components(graph):
        """
        输入k-core子图，统计连通块数量
        Count connected components using BFS
        """
        visited = set()  # 记录访问过的点 Visited nodes
        count = 0        # 连通块计数器 Component counter

        def bfs(start):
            q = deque([start])  # 新建BFS队列 Queue for BFS
            visited.add(start)
            while q:
                u = q.popleft()
                for v in graph[u]:
                    if v not in visited:
                        visited.add(v)
                        q.append(v)

        for u in graph:
            if u not in visited:
                bfs(u)   # 每遇到一个未访问点就开始新连通块 Start a new component
                count += 1

        return count  # 返回连通块数量 Return number of components


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
