#!/usr/bin/env python3
# Auto-generated for 5491311

STUDENT_ID = "5491311"
STUDENT_NAME = "Ruixi Yu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _count_k_cores_in_subgraph(G, nodes_in_subgraph, k):
        if len(nodes_in_subgraph) <= k:
            return 0

        # Step 1: Build the subgraph implicitly by calculating degrees
        subgraph_adj = {}
        subgraph_degrees = {}
        
        # Use a set for efficient neighbor lookup
        nodes_set = set(nodes_in_subgraph)

        for u in nodes_in_subgraph:
            subgraph_adj[u] = []
            degree = 0
            for v_neighbor in G.adj_list[u]:
                if v_neighbor in nodes_set:
                    degree += 1
                    subgraph_adj[u].append(v_neighbor)
            subgraph_degrees[u] = degree
            
        # Step 2: Perform k-core decomposition (peeling)
        q = deque()
        for node in nodes_in_subgraph:
            if subgraph_degrees[node] < k:
                q.append(node)

        removed_nodes = set(q)
        while q:
            u = q.popleft()
            for v_neighbor in subgraph_adj[u]:
                if v_neighbor not in removed_nodes:
                    subgraph_degrees[v_neighbor] -= 1
                    if subgraph_degrees[v_neighbor] < k:
                        q.append(v_neighbor)
                        removed_nodes.add(v_neighbor)

        # The remaining nodes form the k-cores
        core_nodes = nodes_set - removed_nodes
        
        if not core_nodes:
            return 0

        # Step 3: Count connected components in the remaining graph
        visited = set()
        component_count = 0
        for node in core_nodes:
            if node not in visited:
                component_count += 1
                # BFS to find all nodes in the current component
                bfs_q = deque([node])
                visited.add(node)
                while bfs_q:
                    curr = bfs_q.popleft()
                    for neighbor in subgraph_adj[curr]:
                        # Traverse only within the remaining core nodes
                        if neighbor in core_nodes and neighbor not in visited:
                            visited.add(neighbor)
                            bfs_q.append(neighbor)
        
        return component_count


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
        n = G.vertex_num
        tau = [0] * n

        # Iterate through each vertex in the graph G
        for v in range(n):
            # Get the neighbors of v to define the induced subgraph
            neighbors_of_v = G.adj_list[v]

            if len(neighbors_of_v) < k + 1:
                tau[v] = 0
                continue
            
            # Calculate the number of k-cores in the neighbor-induced subgraph
            tau[v] = kCoreBaseStructuralDiversity._count_k_cores_in_subgraph(G, neighbors_of_v, k)
            
        return tau

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
