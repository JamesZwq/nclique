#!/usr/bin/env python3
# Auto-generated for 5509475

STUDENT_ID = "5509475"
STUDENT_NAME = "Danni Cai"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):  #计算每个节点的k-core数量
        total_nodes = G.vertex_num
        diversity_list = [0 for _ in range(total_nodes)]  #每个节点对应的 τ_k(v) 值

        for center_node in range(total_nodes):  #对每个节点操作
            neighbour_list = G.adj_list[center_node]
            if not neighbour_list:
                continue  #如果没有邻居就跳过

            
            induced_subgraph, degrees = kCoreBaseStructuralDiversity._build_induced_subgraph_with_degrees(G, neighbour_list)  #在构建子图的同时返回度数信息

            #计算子图中剩下的k-core连通块数量
            diversity_list[center_node] = kCoreBaseStructuralDiversity._compute_k_core(induced_subgraph, degrees, k)

        return diversity_list

    @staticmethod
    def _build_induced_subgraph_with_degrees(G, nodes):   #构建诱导子图同时计算节点度数
        subgraph = {}  #子图，以邻居为顶点
        degrees = {}    #存储每个节点在诱导子图中的度数
        
        for u in nodes:  #初始化所有节点的邻接表和度数
            subgraph[u] = set()
            degrees[u] = 0
        
        for u in nodes:  #邻接关系，更新度数
            for v in G.adj_list[u]:
                if v in subgraph:  #添加在邻居集中的边
                    if v not in subgraph[u]:
                        subgraph[u].add(v)
                        subgraph[v].add(u)
                        degrees[u] += 1 #重新更新两端的度数
                        degrees[v] += 1
        
        return subgraph, degrees

    @staticmethod
    def _compute_k_core(graph, degrees, k):  #进行k-core剥离,统计剩余节点构成的连通块数量
        if not graph:
            return 0

        
        deg = {node: degrees[node] for node in graph}  #创建度数的可变副本
        removed = set()  #被剥离
        removal_queue = deque()

        for node in graph:  #添加度数不足k的节点
            if deg[node] < k:
                removal_queue.append(node)
                removed.add(node)

        while removal_queue:  #度数不足的节点
            current = removal_queue.popleft()
            for neighbor in graph[current]:
                if neighbor in removed:
                    continue
                    
                deg[neighbor] -= 1  #更新邻居度数
                
                # 检查是否需要移除
                if deg[neighbor] < k:
                    removed.add(neighbor)
                    removal_queue.append(neighbor)

        # 获取剩余节点
        remaining_nodes = [node for node in graph if node not in removed]
        if not remaining_nodes:
            return 0

        return kCoreBaseStructuralDiversity._count_connected_components(graph, removed, remaining_nodes)  #连通分量数量

    @staticmethod
    def _count_connected_components(graph, removed, candidates):  #对剩下的k-core做 BFS
        visited = set()
        count = 0
        
        removed_set = removed  #集合提高查找效率
        
        for start in candidates:  #剩下未删除的节点集合
            if start in visited:
                continue

            queue = deque([start]) #从一个未访问节点开始BFS
            visited.add(start)
            count += 1

            while queue:
                node = queue.popleft()
                for neighbor in graph[node]:
                    if neighbor not in visited and neighbor not in removed_set:
                        visited.add(neighbor)
                        queue.append(neighbor)

        return count

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
