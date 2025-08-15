#!/usr/bin/env python3
# Auto-generated for 5581423

STUDENT_ID = "5581423"
STUDENT_NAME = "Yankai Hou"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(graph_instance, threshold_k):
        total_vertices = graph_instance.vertex_num
        structural_diversity_values = [0] * total_vertices
        
        # 为每个顶点计算结构多样性
        for vertex_id in range(total_vertices):
            adjacent_vertices = set(graph_instance.adj_list[vertex_id])
            
            # 构建邻居顶点的映射关系
            neighbor_vertices = list(adjacent_vertices)
            vertex_index_mapping = {vertex: idx for idx, vertex in enumerate(neighbor_vertices)}
            
            # 构建邻居子图的邻接表
            neighbor_subgraph = [[] for _ in range(len(neighbor_vertices))]
            for idx, current_vertex in enumerate(neighbor_vertices):
                for connected_vertex in graph_instance.adj_list[current_vertex]:
                    if connected_vertex in adjacent_vertices:
                        neighbor_subgraph[idx].append(vertex_index_mapping[connected_vertex])
            
            # 查找k-core组件
            core_components = kCoreBaseStructuralDiversity._extract_k_core_components(neighbor_subgraph, threshold_k)
            structural_diversity_values[vertex_id] = len(core_components)
        
        return structural_diversity_values
    
    @staticmethod
    def _extract_k_core_components(adjacency_structure, k_threshold):
        vertex_count = len(adjacency_structure)
        if vertex_count == 0:
            return []
        
        # 计算核心数
        vertex_core_values = kCoreBaseStructuralDiversity._calculate_vertex_core_numbers(adjacency_structure)
        
        # 筛选满足k-core条件的顶点
        qualified_vertices = [vertex for vertex in range(vertex_count) if vertex_core_values[vertex] >= k_threshold]
        
        if not qualified_vertices:
            return []
        
        # 重新映射顶点索引
        index_remapping = {vertex: new_idx for new_idx, vertex in enumerate(qualified_vertices)}
        remapped_adjacency = [[] for _ in range(len(qualified_vertices))]
        
        for new_idx, original_vertex in enumerate(qualified_vertices):
            for neighbor in adjacency_structure[original_vertex]:
                if neighbor in index_remapping:
                    remapped_adjacency[new_idx].append(index_remapping[neighbor])
        
        # 使用DFS查找连通组件
        exploration_status = [False] * len(qualified_vertices)
        connected_components = []
        
        for start_idx in range(len(qualified_vertices)):
            if not exploration_status[start_idx]:
                current_component = []
                exploration_stack = [start_idx]
                while exploration_stack:
                    current_node = exploration_stack.pop()
                    if not exploration_status[current_node]:
                        exploration_status[current_node] = True
                        current_component.append(qualified_vertices[current_node])
                        for adjacent_node in remapped_adjacency[current_node]:
                            if not exploration_status[adjacent_node]:
                                exploration_stack.append(adjacent_node)
                
                if current_component:
                    connected_components.append(set(current_component))
        
        return connected_components
    
    @staticmethod
    def _calculate_vertex_core_numbers(graph_adjacency):
        node_count = len(graph_adjacency)
        vertex_degrees = [len(graph_adjacency[node]) for node in range(node_count)]
        core_number_array = [0] * node_count
        
        # 找到最大度数
        maximum_degree = max(vertex_degrees) if vertex_degrees else 0
        degree_buckets = [[] for _ in range(maximum_degree + 1)]
        
        # 将顶点按度数分桶
        for node in range(node_count):
            degree_buckets[vertex_degrees[node]].append(node)
        
        # 按度数递增顺序处理顶点
        for current_degree in range(maximum_degree + 1):
            while degree_buckets[current_degree]:
                processed_vertex = degree_buckets[current_degree].pop()
                core_number_array[processed_vertex] = current_degree
                
                # 更新邻居的度数
                for neighbor_vertex in graph_adjacency[processed_vertex]:
                    if vertex_degrees[neighbor_vertex] > current_degree:
                        degree_buckets[vertex_degrees[neighbor_vertex]].remove(neighbor_vertex)
                        vertex_degrees[neighbor_vertex] -= 1
                        degree_buckets[vertex_degrees[neighbor_vertex]].append(neighbor_vertex)
        
        return core_number_array

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
