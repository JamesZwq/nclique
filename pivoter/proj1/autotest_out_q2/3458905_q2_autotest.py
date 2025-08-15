#!/usr/bin/env python3
# Auto-generated for 3458905

STUDENT_ID = "3458905"
STUDENT_NAME = "Tao Yu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules
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
        n = G.vertex_num
        sd = [0] * n

        # if no nodes
        if n == 0:
            return sd

        total_deg=0

        # graph with no edges
        for v in range(n):
          total_deg += len(G.adj_list[v])
        if total_deg == 0:
            return sd

        # maybe k is too small or negative
        if k <= 0:
            return sd

        # calculate core numbers for all nodes
        core_numbers = kCoreBaseStructuralDiversity._calculate_core_numbers(G)

        # find nodes that can be in k-cores
        k_core_vertices = set()
        for v in range(n):
            if core_numbers[v] >= k:
                k_core_vertices.add(v)

        # no k-cores exist
        if not k_core_vertices:
            return sd

        # find connected parts in k-core nodes
        components = kCoreBaseStructuralDiversity._find_connected_components(G, k_core_vertices)

        # process nodes in batches for efficiency
        processed = set()

        # process nodes related to components first
        for component in components:
            kCoreBaseStructuralDiversity._process_component_batch(
                G, component, k, core_numbers, sd, processed
            )

        # make sure we didn't miss any nodes
        for v in range(n):
            if v not in processed:
                sd[v] = kCoreBaseStructuralDiversity._compute_single_vertex_diversity(G, v, k)

        return sd

    @staticmethod
    def _calculate_core_numbers(G):
        """calculate core numbers for all nodes"""
        n = G.vertex_num

        deg = []
        for v in range(n):
            deg.append(len(G.adj_list[v]))

        max_deg = 0
        if deg:
            max_deg = max(deg)
        else:
            max_deg = 0

        # put nodes into buckets by degree
        bins = []
        for i in range(max_deg + 1):
            bins.append([])
        for v in range(n):
            bins[deg[v]].append(v)

        core_numbers = [0] * n
        removed = [False] * n

        # process nodes from low degree to high degree
        for i in range(max_deg + 1):
            while bins[i]:
                v = bins[i].pop()
                if removed[v]:
                    continue

                core_numbers[v] = deg[v]
                removed[v] = True

                # update degrees of neighbors
                for neighbor in G.adj_list[v]:
                    if not removed[neighbor]:
                        old_degree = deg[neighbor]
                        deg[neighbor] -= 1
                        new_degree = deg[neighbor]
                        if new_degree >= 0:
                            bins[new_degree].append(neighbor)

        return core_numbers

    @staticmethod
    def _find_connected_components(G, k_core_nodes):
        """find connected components using BFS"""
        all_components = []
        nodes_visited = set()

        # go through each k-core node
        for start_node in k_core_nodes:
            if start_node in nodes_visited:
                continue

            current_component = []
            search_queue = deque([start_node])
            nodes_visited.add(start_node)

            # BFS to find all nodes in this component
            while len(search_queue) > 0:
                current_node = search_queue.popleft()
                current_component.append(current_node)

                # check all neighbors of current node
                for neighbor in G.adj_list[current_node]:
                    if neighbor in k_core_nodes:
                        if neighbor not in nodes_visited:
                            nodes_visited.add(neighbor)
                            search_queue.append(neighbor)

            all_components.append(current_component)

        return all_components

    @staticmethod
    def _process_component_batch(G, component, k, core_numbers, sd, processed):
        """process nodes related to a component"""

        # collect nodes that need to be processed
        # includes component nodes and their neighbors
        vertices_to_process = set(component)
        for v in component:
            for neighbor in G.adj_list[v]:
                vertices_to_process.add(neighbor)

        # process each node
        for v in vertices_to_process:
            if v in processed:
                continue

            processed.add(v)

            neighbors = G.adj_list[v]

            # isolated node
            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # count neighbors with high core numbers
            high_core_neighbors = []
            for neighbor in neighbors:
                if core_numbers[neighbor] >= k:
                    high_core_neighbors.append(neighbor)

            # not enough high-core neighbors
            if len(high_core_neighbors) < k:
                sd[v] = 0
                continue

            # compute k-cores in neighbor subgraph
            sd[v] = kCoreBaseStructuralDiversity._compute_k_core_count(G, neighbors, k)

    @staticmethod
    def _compute_single_vertex_diversity(G, v, k):
        """compute diversity for one node"""
        neighbors = G.adj_list[v]

        # no neighbors
        if len(neighbors) == 0:
            return 0

        # not enough neighbors
        if len(neighbors) < k:
            return 0

        return kCoreBaseStructuralDiversity._compute_k_core_count(G, neighbors, k)

    @staticmethod
    def _compute_k_core_count(G, neighbors, k):
        """count k-cores in neighbor-induced subgraph"""
        if len(neighbors) < k:
            return 0

        # create mapping from original node id to subgraph index
        vertex_to_idx = {}
        for i in range(len(neighbors)):
            vertex_to_idx[neighbors[i]] = i

        n = len(neighbors)
        adj_list = []
        for i in range(n):
            adj_list.append([])

        # build adjacency list for subgraph
        for i in range(len(neighbors)):
            u = neighbors[i]
            for neighbor_of_u in G.adj_list[u]:
                if neighbor_of_u in vertex_to_idx:
                    j = vertex_to_idx[neighbor_of_u]
                    adj_list[i].append(j)

        # remove nodes with degree < k iteratively
        remaining = set()
        for v in range(n):
            remaining.add(v)

        changed = True
        while changed:
            changed = False
            to_remove = []

            for v in remaining:
                # count degree in remaining subgraph
                current_degree = 0
                for neighbor in adj_list[v]:
                    if neighbor in remaining:
                        current_degree += 1

                if current_degree < k:
                    to_remove.append(v)

            # remove nodes with low degree
            if len(to_remove) > 0:
                changed = True
                for node_to_delete in to_remove:
                    remaining.remove(node_to_delete)

        # check if any nodes left for k-core
        if len(remaining) == 0:
            return 0

        # find how many separate groups exist
        node_visited = set()
        component_count = 0

        # check each remaining node
        for node in remaining:
            if node in node_visited:
                continue

            # start BFS from this node
            bfs_queue = deque([node])
            node_visited.add(node)

            # explore all connected nodes
            while len(bfs_queue) > 0:
                current_node = bfs_queue.popleft()
                for next_node in adj_list[current_node]:
                    if next_node in remaining:
                        if next_node not in node_visited:
                            node_visited.add(next_node)
                            bfs_queue.append(next_node)

            component_count = component_count + 1

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
