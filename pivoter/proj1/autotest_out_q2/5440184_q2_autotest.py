#!/usr/bin/env python3
# Auto-generated for 5440184

STUDENT_ID = "5440184"
STUDENT_NAME = "Yuanbo Li"

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

        # 最终的结果列表，初始化为0
        sd = [0] * G.vertex_num

        # 对图中的每一个顶点 v，计算其 τ_k(v)
        for v in range(G.vertex_num):
            sd[v] = kCoreBaseStructuralDiversity._compute_tau_for_vertex(G, v, k)
        
        return sd

    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def _compute_tau_for_vertex(G, v, k):

 
        neighbors = set(G.adj_list[v])

     
        if len(neighbors) < k:
            return 0
        
        if k == 0:
            return kCoreBaseStructuralDiversity._count_connected_components(G, neighbors)

       
        subgraph_degrees = {}
        for u in neighbors:
            degree = 0
            for neighbor_of_u in G.adj_list[u]:
                if neighbor_of_u in neighbors:
                    degree += 1
            subgraph_degrees[u] = degree

        to_peel = deque([u for u, deg in subgraph_degrees.items() if deg < k])
        
       
        peeled_vertices = set(to_peel)

     
        while to_peel:
            u = to_peel.popleft()
            
        
            for neighbor_of_u in G.adj_list[u]:
         
                if neighbor_of_u in subgraph_degrees and neighbor_of_u not in peeled_vertices:
                    subgraph_degrees[neighbor_of_u] -= 1
                  
                    if subgraph_degrees[neighbor_of_u] < k:
                        peeled_vertices.add(neighbor_of_u)
                        to_peel.append(neighbor_of_u)
        
   
        survivors = neighbors - peeled_vertices
        return kCoreBaseStructuralDiversity._count_connected_components(G, survivors)

    @staticmethod
    def _count_connected_components(G, nodes_to_check):

        if not nodes_to_check:
            return 0

        count = 0
        visited = set()
        
 
        for node in nodes_to_check:
            if node not in visited:
                count += 1
                q = deque([node])
                visited.add(node)
                while q:
                    u = q.popleft()
                 
                    for v_neighbor in G.adj_list[u]:
                    
                        if v_neighbor in nodes_to_check and v_neighbor not in visited:
                            visited.add(v_neighbor)
                            q.append(v_neighbor)
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
