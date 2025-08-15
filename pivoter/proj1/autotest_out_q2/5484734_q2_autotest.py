#!/usr/bin/env python3
# Auto-generated for 5484734

STUDENT_ID = "5484734"
STUDENT_NAME = "Yude Qian"

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
        n = G.vertex_num
        sd = [0] * n

        for v in range(n):
            neighbors = set(G.adj_list[v])

            if len(neighbors) == 0:
                sd[v] = 0
                continue

            neighbor_list = list(neighbors)
            vertex_to_local = {vertex: i for i, vertex in enumerate(neighbor_list)}
            local_adj = [[] for _ in range(len(neighbor_list))]

            for i, u in enumerate(neighbor_list):
                for neighbor in G.adj_list[u]:
                    if neighbor in neighbors:
                        local_adj[i].append(vertex_to_local[neighbor])

            k_cores = kCoreBaseStructuralDiversity._find_k_cores(local_adj, k)
            sd[v] = len(k_cores)

        return sd

    @staticmethod
    def _find_k_cores(adj_list, k):
        """
        Find all k-cores in a graph using iterative vertex removal
        Returns a list of k-cores, where each k-core is a set of vertices
        """
        n = len(adj_list)
        if n == 0:
            return []

        sub_G = [list(neighbors) for neighbors in adj_list]
        remaining = set(range(n))
        k_cores = []

        while remaining:
            deleted = set()
            changed = True

            while changed:
                changed = False
                vertices_to_remove = []

                for u in remaining:
                    if u not in deleted:
                        current_degree = sum(1 for neighbor in sub_G[u]
                                           if neighbor in remaining and neighbor not in deleted)
                        if current_degree < k:
                            vertices_to_remove.append(u)
                            changed = True

                for u in vertices_to_remove:
                    deleted.add(u)

            k_core_vertices = [u for u in remaining if u not in deleted]

            if not k_core_vertices:
                break

            visited = set()

            for start_node in k_core_vertices:
                if start_node not in visited:
                    # BFS to find connected component
                    component = set()
                    queue = deque([start_node])
                    visited.add(start_node)

                    while queue:
                        u = queue.popleft()
                        component.add(u)

                        for neighbor in sub_G[u]:
                            if (neighbor in k_core_vertices and
                                neighbor not in visited and
                                neighbor not in deleted):
                                visited.add(neighbor)
                                queue.append(neighbor)

                    if component:
                        k_cores.append(component)

            remaining = remaining - set(k_core_vertices) - deleted

        return k_cores

    @staticmethod
    def _count_k_cores_alternative(adj_list, k):
        n = len(adj_list)
        if n == 0:
            return 0

        remaining = set(range(n))
        k_core_count = 0

        while remaining:
            sub_G = {}
            for u in remaining:
                sub_G[u] = [v for v in adj_list[u] if v in remaining]

            deleted = set()
            vertices = deque(remaining)

            while vertices:
                u = vertices.popleft()
                if u in deleted:
                    continue

                current_neighbors = [v for v in sub_G[u] if v not in deleted]

                if len(current_neighbors) < k:
                    deleted.add(u)
                    for v in current_neighbors:
                        if v not in deleted:
                            vertices.append(v)

            k_core_vertices = remaining - deleted

            if not k_core_vertices:
                break

            visited = set()

            for node in k_core_vertices:
                if node not in visited:
                    k_core_count += 1
                    # BFS to mark all vertices in this component
                    queue = deque([node])
                    visited.add(node)

                    while queue:
                        current = queue.popleft()
                        for neighbor in sub_G[current]:
                            if neighbor in k_core_vertices and neighbor not in visited:
                                visited.add(neighbor)
                                queue.append(neighbor)

            remaining = remaining - k_core_vertices

        return k_core_count

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
