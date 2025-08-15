#!/usr/bin/env python3
# Auto-generated for 5421856

STUDENT_ID = "5421856"
STUDENT_NAME = "Terrance Ho"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def build_subgraph_neighbour(G, neighbours):
        """
            Build subgraph induced by neighbours
            Given a graph and neighbours
            Return subgraph with the neighbours and edges between them

            time: O(dev(v)^2) -
        """
        subgraph = {}
        neighbour_list = list(neighbours)

        for node in neighbour_list:
            subgraph[node] = []

        # For every edge, add edges that connect vertices within the neighbours set
        for node in neighbour_list:
            for neighbour in G.adj_list[node]:
                if neighbour in neighbours:
                    subgraph[node].append(neighbour)

        return subgraph

    @staticmethod
    def find_connected_components(G):
        """
            Find all connected components in the graph
            Given a graph,
            Return all the different connected components within it

            O(deg(v)^2)
        """

        if not G:
            return []

        visited = set()
        components = []

        for s in G:
            if s not in visited:
                component = []
                queue = deque([s])
                visited.add(s)

                while queue:
                    node = queue.pop()
                    component.append(node)

                    for neighbour in G[node]:
                        if neighbour not in visited:
                            visited.add(neighbour)
                            queue.append(neighbour)
                if component:
                    components.append(component)

        return components

    @staticmethod
    def count_k_cores(subgraph, neighbours, k):
        """
            Count the number of k-cores in the subgraph
            Given a subgraph, neighbours and k
            Return the number of k-cores

            time: O(deg(v)^3)
        """
        if len(neighbours) < k:
            return 0

        nodes_left = set(neighbours)
        k_cores = []

        while len(nodes_left) >= k:

            current_subgraph = {
                node: set(n for n in subgraph[node] if n in nodes_left)
                for node in nodes_left
            }

            changed = True
            while changed:
                changed = False
                remove = [node for node, nbrs in current_subgraph.items() if len(nbrs) < k]

                if remove:
                    changed = True
                    for node in remove:
                        for neighbour in current_subgraph[node]:
                            if neighbour in current_subgraph:
                                current_subgraph[neighbour].remove(node)
                        del current_subgraph[node]

            # if non-empty k-core
            if current_subgraph:
                components = kCoreBaseStructuralDiversity.find_connected_components(current_subgraph)
                k_cores.extend(components)

                # remove all nodes in found k-cores from nodes left
                for component in components:
                    nodes_left -= set(component)
            else:
                break
        return len(k_cores)

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

        # For each vertex v
        for v in range(n):

            neighbours = set(G.adj_list[v])

            # If v has fewer than k neighbors, no k-cores possible
            if len(neighbors) < k:
                tau[v] = 0
                continue

            # Build neighbor-induced subgraph
            neighbor_subgraph = kCoreBaseStructuralDiversity.build_subgraph_neighbour(G, neighbours)

            # Find all k-cores in the neighbor subgraph
            tau[v] = kCoreBaseStructuralDiversity.count_k_cores(neighbor_subgraph, neighbours, k)

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
