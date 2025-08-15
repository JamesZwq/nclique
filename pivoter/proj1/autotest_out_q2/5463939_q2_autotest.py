#!/usr/bin/env python3
# Auto-generated for 5463939

STUDENT_ID = "5463939"
STUDENT_NAME = "Anli Sha"

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
        for k-core based strutural diversity for all vertices
        1. find neighborhood N(v)
        2. neighbor-induced subgraph G[N(v)]
        3. find all k-cores in this subgraph
        4. count the number of k-cores(thy must be connected components)
        """
        
        n = G.vertex_num
        sd = [0] * n
        
        for v in range(n):
            # step1 G_adj_list[5] = [0, 2, 3, 4, 6, 7] -> set
            neighbors = set(G.adj_list[v])
            if len(neighbors) == 0: 
                sd[v] = 0
                continue
            
            # step2 neighbor-induced subgraph
            neighbor_list = list(neighbors) # transform it to a list to make it
            neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbor_list)}
            
            # build adj list for neighbor-induced subgraph
            nbr_adj = [[] for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list): # i is the idx, and u is the neighbor_node
                for w in G.adj_list[u]:
                    if w in neighbor_to_idx and w != u:
                        j = neighbor_to_idx[w]
                        nbr_adj[i].append(j) # here i and j are actually the mapping idx
                        
            sd[v] = kCoreBaseStructuralDiversity._count_k_cores(nbr_adj, k)
        return sd
    
    @staticmethod
    def _count_k_cores(adj_list, k):
        # means all the connected components
        if k == 0:
            return kCoreBaseStructuralDiversity._count_connected_components(adj_list)
        
        n = len(adj_list)
        if n == 0:
            return 0
        
        graph = [set(neighbors) for neighbors in adj_list]
        degrees = [len(neighbors) for neighbors in graph]
        removed = [False] * n
        queue = deque()
        
        for i in range(n):
            if degrees[i] < k:
                queue.append(i)
        
        while queue:
            u = queue.popleft()
            if removed[u]:
                continue
            removed[u] = True
            
            #update its neighbors
            for v in list(graph[u]):
                if not removed[v]:
                    # Remove edge u-v
                    graph[v].discard(u)
                    degrees[v] -= 1
                    # if it's degree drops below k, mark for removal
                    if degrees[v] < k:
                        queue.append(v)
            
            graph[u].clear()
            degrees[u] = 0
        
        # Count connected components among remaining vertices (integrated function)
        return kCoreBaseStructuralDiversity._count_connected_components(graph, removed_nodes=removed)
            
    @staticmethod
    def _count_connected_components(adj_list, removed_nodes=None):
        """
        - removed_nodes=None: count all components (k=0 case)
        - removed_nodes provided: skip removed vertices (k>0 case)
        """
        n = len(adj_list)
        visited = [False] * n
        components = 0
        
        for i in range(n):
            # Skip if node is removed or already visited
            if (removed_nodes and removed_nodes[i]) or visited[i]:
                continue
                
            components += 1
            # BFS to mark all components (or "Also use BFS" for remaining components)
            queue = deque([i])
            visited[i] = True
            
            
            while queue:
                u = queue.popleft()
                for v in adj_list[u]:
                    # Skip removed or visited vertices
                    if (removed_nodes and removed_nodes[v]) or visited[v]:
                        continue
                    visited[v] = True
                    queue.append(v)
                            
        return components

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
