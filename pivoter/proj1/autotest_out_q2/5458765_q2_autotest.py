#!/usr/bin/env python3
# Auto-generated for 5458765

STUDENT_ID = "5458765"
STUDENT_NAME = "Yi Fan"

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

        # Iterate through each vertex in the graph
        for idx in range(n):
            # Fetch the neighbours of the current vertex
            neighbours = G.adj_list[idx]
            # Skip if no neighbours
            if not neighbours:
                continue

            # Construct a subgraph with the neighbours
            subgraph = kCoreBaseStructuralDiversity.construct_subgraph(G, idx)

            # Compute the k-core subgraph
            k_core_subgraph = kCoreBaseStructuralDiversity.online_compute_core(subgraph, k)

            # Count the number of components in the k-core subgraph
            count = kCoreBaseStructuralDiversity.components_count(k_core_subgraph)
             
            sd[idx] = count # Assign the value of τ_k(v)
        return sd
    
    @staticmethod
    def construct_subgraph(G, vertex):
        """
        construct induced subgraph of the neighbors of vertex
        """
        neighbours = G.adj_list[vertex]

        # Create a mapping from original vertex IDs to indices in the subgraph
        vertex_idx_mapping = {u:i for i,u in enumerate(neighbours)}

        # Initialize the subgraph
        subgraph = [list() for _ in range(len(neighbours))]

        for u in neighbours:
            u_idx = vertex_idx_mapping[u]
            for v in G.adj_list[u]:
                if v in vertex_idx_mapping:
                    v_idx = vertex_idx_mapping[v]
                    subgraph[u_idx].append(v_idx)
        
        return subgraph
    
    @staticmethod
    def online_compute_core(subgraph,k):
        """ 
        k-core online computation, iteratively removing vertices with degree < k
        """
        length_subgraph = len(subgraph)
        degree = [len(subgraph[i]) for i in range(length_subgraph)]
        core = [0]*length_subgraph # 1 deleted, 0 not deleted

        # Initialize, add all vertices with degree < k to the queue
        queue = deque(i for i in range(length_subgraph) if degree[i] < k)
        
        for i in queue:
            core[i] = 1
        
        # Remove vertices iteratively and update neighbour degrees
        while queue:
            u = queue.popleft()
            for v in subgraph[u]:
                if not core[v]: # Only process unremoved vertices
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)
                        core[v] = 1
        
        # Construct remaining subgraph (only vertices with degree >= k)
        remaing_vertices = [i for i in range(length_subgraph) if not core[i]]
        new_graph=[list() for _ in range(len(remaing_vertices))]

        # Remap old vertex indices to new indices in the reduced graph
        old2new_mapping = {old:new for new, old in enumerate(remaing_vertices)}

        # Rebuild adjecency list of the remaining k-core subgraph
        for old in remaing_vertices:
            new = old2new_mapping[old]
            for neighbour in subgraph[old]:
                if neighbour in old2new_mapping:
                    new_graph[new].append(old2new_mapping[neighbour])
        
        return new_graph

    @staticmethod
    def components_count(subgraph):
        """ 
        Count connected components in an undirected graph
        """
        length_subgraph = len(subgraph)
        visited = set()
        components = 0
        
        # Traverse all vertices in the subgraph BFS
        for i in range(length_subgraph):
            # Only unvisited and has neighbours vertices
            if i not in visited and subgraph[i]:
                components += 1
                queue = deque([i])
                # BFS to mark all reachable vertices from i
                visited.add(i)
                while queue:
                    current = queue.popleft()
                    for neighbour in subgraph[current]:
                        if neighbour not in visited:
                            visited.add(neighbour)
                            queue.append(neighbour)
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
