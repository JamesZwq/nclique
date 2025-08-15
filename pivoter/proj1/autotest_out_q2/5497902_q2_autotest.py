#!/usr/bin/env python3
# Auto-generated for 5497902

STUDENT_ID = "5497902"
STUDENT_NAME = "Silin Zhou"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def _count_components(adj_list, nodes_to_consider):
        """
        Count the number of connected components in the given set of nodes using BFS.
        adj_list: adjacency list of the subgraph.
        nodes_to_consider: find connected components within this set of nodes.
        """
        if not nodes_to_consider:
            return 0

        count = 0
        visited = set()

        for node in nodes_to_consider:
            if node not in visited:
                count += 1
                q = deque([node])
                visited.add(node)
                while q:
                    curr = q.popleft()
                    if curr in adj_list:
                        for neighbor in adj_list[curr]:
                            if neighbor in nodes_to_consider and neighbor not in visited:
                                visited.add(neighbor)
                                q.append(neighbor)
        return count

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
        num_vertices = G.vertex_num
        structural_diversities = [0] * num_vertices

        for v in range(num_vertices):

            neighbor_set = set(G.adj_list[v])

            if len(neighbor_set) < k:
                continue

            subgraph_adj = {}
            for u in neighbor_set:
                u_internal_neighbors = [n for n in G.adj_list[u] if n in neighbor_set]
                subgraph_adj[u] = u_internal_neighbors

            subgraph_degrees = {u: len(adj) for u, adj in subgraph_adj.items()}

            if k == 0:
                structural_diversities[v] = kCoreBaseStructuralDiversity._count_components(subgraph_adj, neighbor_set)
                continue

            q = deque([u for u, deg in subgraph_degrees.items() if deg < k])
            removed_nodes = set(q)

            while q:
                u = q.popleft()

                for neighbor in subgraph_adj[u]:
                    if neighbor not in removed_nodes:
                        subgraph_degrees[neighbor] -= 1
                        if subgraph_degrees[neighbor] < k:
                            removed_nodes.add(neighbor)
                            q.append(neighbor)

            k_core_nodes = neighbor_set - removed_nodes

            num_components = kCoreBaseStructuralDiversity._count_components(subgraph_adj, k_core_nodes)
            structural_diversities[v] = num_components

        return structural_diversities


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
