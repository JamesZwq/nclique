#!/usr/bin/env python3
# Auto-generated for 5485899

STUDENT_ID = "5485899"
STUDENT_NAME = "Jiangzihan Sun"

# ======= 学生代码 =======
from collections import defaultdict, deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        # tau 是每个节点的结构多样性分数
        tau = [0 for _ in range(G.vertex_num)]

        for v in range(G.vertex_num):
            # 获取该节点的邻居集合 N(v)
            neighbors = G.adj_list[v]

            # 创建一个集合来保存诱导子图中的节点
            induced_node_set = set()
            for u in neighbors:
                induced_node_set.add(u)

            # 创建一个邻接表来表示诱导子图
            induced_subgraph = defaultdict(set)

            # 遍历诱导节点集合，构建它们之间的边
            for u in induced_node_set:
                for w in G.adj_list[u]:
                    if w in induced_node_set:
                        induced_subgraph[u].add(w)

            # 调用函数计算该诱导子图中的 k-core 连通分量个数
            k_core_components = kCoreBaseStructuralDiversity._count_k_core_components(induced_subgraph, k)

            tau[v] = k_core_components

        return tau

    @staticmethod
    def _count_k_core_components(graph, k):
        # 如果诱导子图是空的，直接返回0
        if len(graph) == 0:
            return 0

        # 创建字典记录每个节点的度数
        degree_dict = {}
        for node in graph:
            degree_dict[node] = len(graph[node])

        # removed 集合记录被删除的节点
        removed = set()

        # 循环进行度数剥离操作，直到图稳定
        changed = True
        while changed:
            changed = False
            # 遍历图中的每一个节点
            for node in list(graph.keys()):
                if node in removed:
                    continue
                current_degree = degree_dict[node]
                if current_degree < k:
                    # 标记该节点为删除
                    removed.add(node)
                    # 更新邻居节点的度数
                    for neighbor in graph[node]:
                        if neighbor not in removed:
                            degree_dict[neighbor] -= 1
                    changed = True  # 本轮发生了变化，需要继续

        # 统计剩余图中连通块的数量
        visited = set()
        component_count = 0

        # 遍历所有剩下的节点，用 BFS 找连通块
        for node in graph:
            if node in removed:
                continue
            if node in visited:
                continue

            # 使用队列来进行 BFS
            bfs_queue = deque()
            bfs_queue.append(node)
            visited.add(node)

            while len(bfs_queue) > 0:
                current = bfs_queue.popleft()
                for neighbor in graph[current]:
                    if neighbor not in removed and neighbor not in visited:
                        visited.add(neighbor)
                        bfs_queue.append(neighbor)

            # 每完成一次 BFS，就找到一个连通分量
            component_count += 1

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
