#!/usr/bin/env python3
# Auto-generated for 5528352

STUDENT_ID = "5528352"
STUDENT_NAME = "Yili Sun"

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
        # n = G.vertex_num
        # sd = [0] * n
        # return sd
        vertex_count = G.vertex_num
        sd = [0] * vertex_count
        
        for vertex in range(vertex_count):
            adjacent_nodes = G.adj_list[vertex]
            if not adjacent_nodes:
                sd[vertex] = 0
                continue
                
            node_neighborhood = set(adjacent_nodes)
            
            local_graph = {}
            for node in adjacent_nodes:
                local_graph[node] = []
                for connected_node in G.adj_list[node]:
                    if connected_node in node_neighborhood:
                        local_graph[node].append(connected_node)
            
            sd[vertex] = kCoreBaseStructuralDiversity._compute_k_core_components(local_graph, k)
        
        return sd

    @staticmethod
    def _compute_k_core_components(graph_data, core_threshold):
        """
        Calculate the number of k-core components in a graph.
        """
        if not graph_data or core_threshold < 0:
            return 0
        
        if core_threshold == 0:
            return kCoreBaseStructuralDiversity._get_component_count(graph_data)
        
        current_graph = {v: set(edges) for v, edges in graph_data.items()}
        
        modified = True
        while modified:
            modified = False
            removal_candidates = []
            
            for node in current_graph:
                if len(current_graph[node]) < core_threshold:
                    removal_candidates.append(node)
            
            for node in removal_candidates:
                for neighbor in current_graph[node]:
                    current_graph[neighbor].discard(node)
                current_graph.pop(node)
                modified = True
        
        if not current_graph:
            return 0
        
        explored_nodes = set()
        component_count = 0
        
        for node in current_graph:
            if node not in explored_nodes:
                node_queue = deque([node])
                explored_nodes.add(node)
                
                while node_queue:
                    current_node = node_queue.popleft()
                    for adjacent in current_graph[current_node]:
                        if adjacent not in explored_nodes:
                            explored_nodes.add(adjacent)
                            node_queue.append(adjacent)
                
                component_count += 1
        
        return component_count
    
    @staticmethod
    def _get_component_count(graph_structure):
        """
        Determine the number of connected components in a graph.
        """
        if not graph_structure:
            return 0
        
        marked_nodes = set()
        total_components = 0
        
        for node in graph_structure:
            if node not in marked_nodes:
                traversal_queue = deque([node])
                marked_nodes.add(node)
                
                while traversal_queue:
                    current_vertex = traversal_queue.popleft()
                    for connected_vertex in graph_structure[current_vertex]:
                        if connected_vertex not in marked_nodes:
                            marked_nodes.add(connected_vertex)
                            traversal_queue.append(connected_vertex)
                
                total_components += 1
        
        return total_components


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
