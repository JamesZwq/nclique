#!/usr/bin/env python3
# Auto-generated for 5415372

STUDENT_ID = "5415372"
STUDENT_NAME = "Ruoyao Xu"

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
        # Calculate structural diversity for each vertex
        for vertex in range(n):
            neighbors = G.adj_list[vertex]

            # If no neighbors, diversity is 0
            if len(neighbors) == 0:
                sd[vertex] = 0
                continue

            # Build neighbor-induced subgraph
            neighbor_subgraph = kCoreBaseStructuralDiversity.build_neighbor_subgraph(G, vertex)


            # Count number of k-cores in neighbor subgraph
            k_core_count = kCoreBaseStructuralDiversity.count_k_cores_in_subgraph(neighbor_subgraph, k)

            sd[vertex] = k_core_count

        return sd

    @staticmethod
    def build_neighbor_subgraph(G, vertex):

        #Build neighbor-induced subgraph for specified vertex

        #Return adjacency list of neighbor subgraph
        neighbors = G.adj_list[vertex]
        neighbor_count = len(neighbors)


        # If no neighbors, return empty graph
        if neighbor_count == 0:
            return []


        # Create mapping from neighbors to new IDs
        neighbor_to_id = {}
        for i in range(neighbor_count):
            original_neighbor = neighbors[i]
            neighbor_to_id[original_neighbor] = i


        # Build adjacency list of neighbor subgraph
        subgraph_adj = []
        for i in range(neighbor_count):
            subgraph_adj.append([])


        # Add edges between neighbors
        for i in range(neighbor_count):
            current_neighbor = neighbors[i]


            # Check all adjacent vertices of current neighbor
            for adjacent_vertex in G.adj_list[current_neighbor]:

                # If adjacent vertex is also in neighbor list, add edge
                if adjacent_vertex in neighbor_to_id:
                    target_id = neighbor_to_id[adjacent_vertex]
                    subgraph_adj[i].append(target_id)

        return subgraph_adj

    @staticmethod
    def count_k_cores_in_subgraph(subgraph_adj, k):


        #Count number of k-core connected components in subgraph

        subgraph_size = len(subgraph_adj)


        # Empty graph has no k-core
        if subgraph_size == 0:
            return 0


        # Perform k-core decomposition to find vertices in k-core
        k_core_vertices = kCoreBaseStructuralDiversity.find_k_core_vertices(subgraph_adj, k)


        # If no k-core vertices
        if len(k_core_vertices) == 0:
            return 0


        # Count connected components of k-core vertices
        component_count = kCoreBaseStructuralDiversity.count_connected_components(subgraph_adj, k_core_vertices)

        return component_count

    @staticmethod
    def find_k_core_vertices(subgraph_adj, k):

        # Use simple iterative method to find vertices in k-core

        subgraph_size = len(subgraph_adj)


        # Calculate degree of each vertex
        degrees = []
        for i in range(subgraph_size):
            degree = len(subgraph_adj[i])
            degrees.append(degree)


        # Mark removed vertices
        is_removed = [False] * subgraph_size


        # Repeatedly remove vertices with degree < k
        something_changed = True
        while something_changed:
            something_changed = False


            # Check each non-removed vertex
            for vertex in range(subgraph_size):
                if is_removed[vertex]:
                    continue


                # If degree < k, remove this vertex
                if degrees[vertex] < k:
                    is_removed[vertex] = True
                    something_changed = True


                    # Update degrees of all neighbors
                    for neighbor in subgraph_adj[vertex]:
                        if not is_removed[neighbor]:
                            degrees[neighbor] = degrees[neighbor] - 1


        # Collect remaining vertices (vertices in k-core)
        k_core_vertices = []
        for vertex in range(subgraph_size):
            if not is_removed[vertex]:
                k_core_vertices.append(vertex)

        return k_core_vertices

    @staticmethod
    def count_connected_components(subgraph_adj, target_vertices):

        #Count number of connected components in specified vertex set

        if len(target_vertices) == 0:
            return 0


        # Convert target vertices to set for fast lookup
        target_set = set()
        for vertex in target_vertices:
            target_set.add(vertex)


        # Record visited vertices
        visited = set()
        component_count = 0


        # Start BFS for each unvisited target vertex
        for start_vertex in target_vertices:
            if start_vertex in visited:
                continue


            # Found a new connected component
            component_count = component_count + 1


            # BFS traversal of current connected component
            queue = deque()
            queue.append(start_vertex)
            visited.add(start_vertex)

            while len(queue) > 0:
                current_vertex = queue.popleft()


                # Check all neighbors of current vertex
                for neighbor in subgraph_adj[current_vertex]:

                    # If neighbor is in target set and not visited
                    if neighbor in target_set and neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)

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
