#!/usr/bin/env python3
# Auto-generated for 5435628

STUDENT_ID = "5435628"
STUDENT_NAME = "Shuo Qiu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity:
    """
    Implements the k-core based structural diversity algorithm.
    """
    @staticmethod
    def process(G, k):
        """
        Computes the k-core based structural diversity τ_k(v) for every vertex v.
        
        Time Complexity: O(V * (avg_deg + avg_deg_neighbors)) which simplifies
                         to O(V*E) in worst case, but much faster on real graphs.
        - Outer loop: V iterations.
        - Inside loop:
            - Building subgraph: O(deg(v) + sum(deg(u))) for u in N(v).
            - k-core decomposition: Proportional to edges in subgraph.
            - CC counting: Proportional to nodes+edges in k-core subgraph.
        """
        tau_values = [0] * G.vertex_num
        for v in range(G.vertex_num):
            tau_values[v] = kCoreBaseStructuralDiversity._compute_diversity_for_vertex(G, v, k)
        return tau_values

    @staticmethod
    def _compute_diversity_for_vertex(G, v_id, k):
        """Computes τ_k for a single vertex."""
        # Step 1: Get the set of neighbors for the vertex v_id.
        neighbors = set(G.adj_list[v_id])
        if not neighbors:
            return 0
            
        # Step 2: Build the neighbor-induced subgraph.
        subgraph = kCoreBaseStructuralDiversity._create_induced_subgraph(G, neighbors)

        # Step 3: Find all vertices belonging to the k-core of the subgraph.
        k_core_nodes = kCoreBaseStructuralDiversity._get_k_core_nodes(subgraph, k)

        if not k_core_nodes:
            return 0

        # Step 4: Count the connected components within the k-core subgraph.
        # We can reuse the original subgraph structure and just provide the valid nodes.
        return kCoreBaseStructuralDiversity._count_connected_components(subgraph, k_core_nodes)

    @staticmethod
    def _create_induced_subgraph(G, node_set):
        """Creates a subgraph induced by a given set of nodes."""
        subgraph = defaultdict(list)
        for u in node_set:
            for v_neighbor in G.adj_list[u]:
                if v_neighbor in node_set:
                    subgraph[u].append(v_neighbor)
        return subgraph
    
    @staticmethod
    def _get_k_core_nodes(subgraph, k):
        """
        Finds the set of nodes that form the k-core of a given subgraph.
        This is implemented using the standard iterative peeling algorithm.
        Fixed version that correctly handles the peeling process.
        """
        if not subgraph:
            return set()
        
        # Create a working copy of the graph
        remaining_nodes = set(subgraph.keys())
        degrees = {u: len(neighbors) for u, neighbors in subgraph.items()}
        
        # Keep removing nodes with degree < k until no more can be removed
        changed = True
        while changed:
            changed = False
            to_remove = []
            
            # Find all nodes with degree < k
            for node in remaining_nodes:
                current_degree = 0
                # Count degree among remaining nodes only
                for neighbor in subgraph[node]:
                    if neighbor in remaining_nodes:
                        current_degree += 1
                
                if current_degree < k:
                    to_remove.append(node)
            
            # Remove the nodes
            if to_remove:
                changed = True
                for node in to_remove:
                    remaining_nodes.discard(node)
        
        return remaining_nodes

    @staticmethod
    def _count_connected_components(graph, nodes_to_check):
        """Counts the number of connected components in a graph on a given subset of its nodes."""
        if not nodes_to_check:
            return 0
        
        visited = set()
        count = 0
        
        for node in nodes_to_check:
            if node not in visited:
                count += 1
                q = deque([node])
                visited.add(node)
                while q:
                    u = q.popleft()
                    # Iterate through neighbors in the original graph structure
                    for v_neighbor in graph.get(u, []):
                        # Only consider neighbors that are part of the valid node set
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
