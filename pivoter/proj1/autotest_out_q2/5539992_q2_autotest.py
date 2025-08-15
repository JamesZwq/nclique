#!/usr/bin/env python3
# Auto-generated for 5539992

STUDENT_ID = "5539992"
STUDENT_NAME = "Qianhe Guo"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(graph, k_val):
        """
        Compute k-core structural diversity τₖ(v) for each node.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        Returns
        -------
        List[int]  # τ_k(v) for all v
        """
        num_nodes = graph.vertex_num
        diversity_result = [0] * num_nodes  # Store τₖ(v) for each node

        for node_id in range(num_nodes):
            neighbors = graph.adj_list[node_id]
            if not neighbors:
                diversity_result[node_id] = 0
                continue

            # Build neighbor-induced subgraph (excluding the node itself)
            unique_neighbors = set(neighbors)
            subgraph = kCoreBaseStructuralDiversity._build_induced_subgraph(unique_neighbors, graph)

            # Count the number of connected k-core components in the subgraph
            diversity_result[node_id] = kCoreBaseStructuralDiversity._count_kcore_components(subgraph, k_val)

        return diversity_result

    @staticmethod
    def _build_induced_subgraph(nodes, original_graph):
        """
        Construct the subgraph induced by a node's neighbors.
        """
        sub_adj_list = {u: [] for u in nodes}
        for u in nodes:
            for v in original_graph.adj_list[u]:
                if v != u and v in nodes:  # Exclude self-loops
                    sub_adj_list[u].append(v)
        return sub_adj_list

    @staticmethod
    def _count_kcore_components(subgraph, k_threshold):
        """
        Prune nodes with degree < k and count connected components in the remaining k-core.
        """
        if not subgraph:
            return 0

        # Step 1: Initialize degrees
        degree = {node: len(adj) for node, adj in subgraph.items()}

        # Step 2: Iteratively remove nodes with degree < k
        removed_nodes = set()
        queue = deque([node for node, deg in degree.items() if deg < k_threshold])
        removed_nodes.update(queue)
        while queue:
            current = queue.popleft()
            for neighbor in subgraph[current]:
                if neighbor not in removed_nodes:
                    degree[neighbor] -= 1
                    if degree[neighbor] < k_threshold:
                        removed_nodes.add(neighbor)
                        queue.append(neighbor)

        # Step 3: Remaining nodes form the k-core
        remaining_nodes = [node for node in subgraph if node not in removed_nodes]
        if not remaining_nodes:
            return 0

        # Step 4: Count connected components in the remaining k-core
        return kCoreBaseStructuralDiversity._count_connected_components(subgraph, remaining_nodes, removed_nodes)

    @staticmethod
    def _count_connected_components(graph, valid_nodes, excluded_nodes):
        """
        BFS to count connected components in the valid k-core nodes.
        """
        visited = set()
        component_count = 0

        for start_node in valid_nodes:
            if start_node not in visited:
                component_count += 1
                queue = deque([start_node])
                visited.add(start_node)
                while queue:
                    current = queue.popleft()
                    for neighbor in graph[current]:
                        if neighbor not in visited and neighbor not in excluded_nodes:
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
