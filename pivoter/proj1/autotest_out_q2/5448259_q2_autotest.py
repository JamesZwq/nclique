#!/usr/bin/env python3
# Auto-generated for 5448259

STUDENT_ID = "5448259"
STUDENT_NAME = "Boheng Hu"

# ======= 学生代码 =======
###############################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _find_connected_components(G, core_nodes):
        """
        A helper function to find the number of connected components on a given set of nodes (the k-core vertices).
        It performs a BFS traversal on the structure of the original graph G, but restricts the traversal to the core_nodes set.

        Args:
            G (UndirectedUnweightedGraph): The original graph object.
            core_nodes (set): The set of vertices in the k-core.

        Returns:
            int: The number of connected components.
        """
        if not core_nodes:
            return 0
        
        visited = set()
        component_count = 0
        
        # Iterate over all nodes belonging to the k-core
        for node in core_nodes:
            if node not in visited:
                # Found an unvisited node, which means a new connected component is found
                component_count += 1
                q = deque([node])
                visited.add(node)
                while q:
                    curr = q.popleft()
                    # Iterate over the neighbors of the current node
                    for neighbor in G.adj_list[curr]:
                        # Only continue if the neighbor is also in the k-core and has not been visited
                        if neighbor in core_nodes and neighbor not in visited:
                            visited.add(neighbor)
                            q.append(neighbor)
        return component_count

    @staticmethod
    def process(G, k):
        """
        Calculates the k-core-based structural diversity τ_k(v) for all vertices in the graph.

        Args:
            G (UndirectedUnweightedGraph): The input undirected, unweighted graph.
            k (int): The k value for the k-core.

        Returns:
            List[int]: A list where the v-th element is the τ_k(v) value for vertex v.
        """
        n = G.vertex_num
        tau_values = [0] * n
        
        # To speed up lookups and intersection operations, pre-convert the adjacency lists to a list of sets.
        adj_sets = [set(neighbors) for neighbors in G.adj_list]

        # Iterate over each vertex v in the graph.
        for v in range(n):
            neighbors_of_v = adj_sets[v]
            
            # If the number of neighbors of v is less than k, its neighbor-induced subgraph cannot form a k-core.
            if len(neighbors_of_v) < k:
                tau_values[v] = 0
                continue

            # --- Step 1: Efficiently calculate the degrees of v's neighbors in its induced subgraph ---
            # The degree in the subgraph is the number of common neighbors with v.
            subgraph_degrees = {u: len(adj_sets[u].intersection(neighbors_of_v)) for u in neighbors_of_v}

            # --- Step 2: Run the k-core peeling algorithm on the neighbor-induced subgraph ---
            q = deque()
            removed_nodes = set()

            # Initialize the queue with all nodes that have a degree less than k in the subgraph.
            for node, degree in subgraph_degrees.items():
                if degree < k:
                    q.append(node)
                    removed_nodes.add(node)
            
            # Start the peeling process.
            while q:
                u = q.popleft()
                
                # u is removed, so we need to update the degrees of its neighbors in the subgraph.
                # A neighbor w of u must be a neighbor of both u and v.
                for w in adj_sets[u].intersection(neighbors_of_v):
                    if w not in removed_nodes:
                        subgraph_degrees[w] -= 1
                        # If the degree of neighbor w also drops below k, add it to the removal queue.
                        if subgraph_degrees[w] < k:
                            q.append(w)
                            removed_nodes.add(w)

            # --- Step 3: Determine the set of vertices in the k-core ---
            # The neighbor nodes that were not removed are the members of the k-core.
            core_nodes = neighbors_of_v - removed_nodes
            
            # --- Step 4: Calculate the number of connected components on the k-core vertex set ---
            tau_values[v] = kCoreBaseStructuralDiversity._find_connected_components(G, core_nodes)
            
        return tau_values


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
