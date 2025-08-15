#!/usr/bin/env python3
# Auto-generated for 5561320

STUDENT_ID = "5561320"
STUDENT_NAME = "Yidi Wu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import defaultdict, deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    # Compute k-core-based structural diversity for all vertices.
    @staticmethod
    def process(G, k):
        n = G.vertex_num
        r = [0] * n
        
        # Convert adjacency lists to sets 
        adj_sets = [set(neighbors) for neighbors in G.adj_list]
        
        # Process each vertex independently
        for v in range(n):
            neighbors = adj_sets[v]
            
            # Skip if not enough neighbors to form k-core
            if len(neighbors) < k:
                continue
                
            # Build neighbor-induced subgraph
            subgraph = kCoreBaseStructuralDiversity._build_neighbor_subgraph(
                neighbors, adj_sets
            )
            
            # Count k-core components in the subgraph
            r[v] = kCoreBaseStructuralDiversity._count_kcore_components(
                subgraph, k
            )
        
        return r
    
    # Build the subgraph induced by the given neighbors.   
    @staticmethod
    def _build_neighbor_subgraph(neighbors, adj_sets):
        subgraph = {u: set() for u in neighbors}
        
        for u in neighbors:
            for w in adj_sets[u]:
                if w in neighbors:
                    subgraph[u].add(w)
                    subgraph[w].add(u)
        
        return subgraph
    
    # Count the number of k-core components in the given subgraph.
    @staticmethod
    def _count_kcore_components(subgraph, k):
        if k == 0:
            # 0-core is just connected components
            return kCoreBaseStructuralDiversity._count_components(
                subgraph, set(subgraph.keys())
            )
        
        # Perform k-core decomposition using peeling algorithm
        remaining = set(subgraph.keys())
        degrees = {u: len(neighbors) for u, neighbors in subgraph.items()}
        
        # Queue for vertices to be removed
        to_remove = deque([u for u, deg in degrees.items() if deg < k])
        
        while to_remove:
            u = to_remove.popleft()
            if u not in remaining:
                continue
                
            remaining.remove(u)
            
            # Update degrees of neighbors
            for neighbor in subgraph[u]:
                if neighbor in remaining:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] == k - 1:
                        to_remove.append(neighbor)
        
        if not remaining:
            return 0
        
        # Count connected components in the remaining k-core
        return kCoreBaseStructuralDiversity._count_components(subgraph, remaining)
    
    # Count connected components using BFS.
    @staticmethod
    def _count_components(graph, vertices):
        visited = set()
        components = 0
        
        for start in vertices:
            if start in visited:
                continue
                
            # BFS to mark all vertices in this component
            components += 1
            queue = deque([start])
            visited.add(start)
            
            while queue:
                current = queue.popleft()
                for neighbor in graph[current]:
                    if neighbor in vertices and neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
        
        return components

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
