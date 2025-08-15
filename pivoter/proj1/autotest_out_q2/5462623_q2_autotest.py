#!/usr/bin/env python3
# Auto-generated for 5462623

STUDENT_ID = "5462623"
STUDENT_NAME = "Yu Xia"

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

        for v in range(n):
            neighbor_nodes = G.adj_list[v]
            if not neighbor_nodes:
                sd[v] = 0
                continue

            induced_subgraph = kCoreBaseStructuralDiversity._build_induced_subgraph(G, neighbor_nodes)

            core_component_count = kCoreBaseStructuralDiversity._compute_core_components(induced_subgraph, k)
            sd[v] = core_component_count

        return sd
    
    @staticmethod
    def _build_induced_subgraph(graph, nodes):
        node_set = set(nodes)
        subgraph = {u: set() for u in node_set}
        for u in node_set:
            for neighbor in graph.adj_list[u]:
                if neighbor in node_set:
                    subgraph[u].add(neighbor)
        return subgraph

    @staticmethod
    def _compute_core_components(subgraph, k):
        if not subgraph:
            return 0

        degree_map = {node: len(neighbors) for node, neighbors in subgraph.items()}

        to_remove = deque(node for node, deg in degree_map.items() if deg < k)
        removed = set(to_remove)

        while to_remove:
            current = to_remove.popleft()
            for neighbor in subgraph[current]:
                if neighbor not in removed:
                    degree_map[neighbor] -= 1
                    if degree_map[neighbor] < k:
                        removed.add(neighbor)
                        to_remove.append(neighbor)

        remaining_nodes = [node for node in subgraph if node not in removed]

        if not remaining_nodes:
            return 0

        visited = set()
        component_count = 0

        for node in remaining_nodes:
            if node not in visited:
                component_count += 1
                kCoreBaseStructuralDiversity._bfs_traverse(subgraph, node, visited, removed)

        return component_count

    @staticmethod
    def _bfs_traverse(graph, start, visited, excluded):
        queue = deque([start])
        visited.add(start)
        while queue:
            u = queue.popleft()
            for v in graph[u]:
                if v not in excluded and v not in visited:
                    visited.add(v)
                    queue.append(v)


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
