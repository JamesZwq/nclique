#!/usr/bin/env python3
# Auto-generated for 5641886

STUDENT_ID = "5641886"
STUDENT_NAME = "Yucheng Zhu"

# ======= 学生代码 =======
from collections import deque

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
        List[int]  # Number of k-core components for each vertex's neighbor-induced subgraph
        """
        # Initialize result array for tau_k(v)
        num_vertices = G.vertex_num
        result = [0] * num_vertices

        # Loop through each vertex to compute k-core components in its neighbor subgraph
        for vertex in range(num_vertices):
            # Get neighbors of current vertex
            neighbor_set = set(G.adj_list[vertex])  # Remove duplicates
            if not neighbor_set:
                result[vertex] = 0
                continue

            # Build subgraph with only neighbors (not including vertex itself)
            subgraph = {u: [] for u in neighbor_set}
            for u in neighbor_set:
                for neighbor in G.adj_list[u]:
                    if neighbor in neighbor_set:
                        subgraph[u].append(neighbor)

            # Compute number of k-core components in subgraph
            result[vertex] = kCoreBaseStructuralDiversity._count_k_core_components(subgraph, k)

        return result
    
    @staticmethod
    def _count_k_core_components(subgraph, k):
        # Handle empty subgraph
        if not subgraph:
            return 0
        
        # Calculate initial degrees for each node in subgraph
        node_degrees = {node: len(neighbors) for node, neighbors in subgraph.items()}
        removed_nodes = set()
        node_queue = deque()

        # Find nodes with degree < k (they can't be in k-core)
        for node in subgraph:
            if node_degrees[node] < k:
                node_queue.append(node)
                removed_nodes.add(node)

        # Remove nodes with degree < k and update neighbors' degrees
        while node_queue:
            current_node = node_queue.popleft()
            for neighbor in subgraph[current_node]:
                if neighbor not in removed_nodes:
                    node_degrees[neighbor] -= 1
                    if node_degrees[neighbor] < k and neighbor not in removed_nodes:
                        node_queue.append(neighbor)
                        removed_nodes.add(neighbor)

        # Get remaining nodes (k-core nodes)
        remaining_nodes = [node for node in subgraph if node not in removed_nodes]
        if not remaining_nodes:
            return 0

        # Count connected components in k-core using BFS
        visited = set()
        component_count = 0
        bfs_queue = deque()

        for node in remaining_nodes:
            if node not in visited:
                component_count += 1
                bfs_queue.append(node)
                visited.add(node)
                while bfs_queue:
                    current = bfs_queue.popleft()
                    for neighbor in subgraph[current]:
                        if neighbor not in removed_nodes and neighbor not in visited:
                            visited.add(neighbor)
                            bfs_queue.append(neighbor)

        return component_count

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
