#!/usr/bin/env python3
# Auto-generated for 5495395

STUDENT_ID = "5495395"
STUDENT_NAME = "Shumin Chen"

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
        tau = [0] * n

        adj_sets = [set(neighbors) for neighbors in G.adj_list]

        for v in range(n):
            neighbors = G.adj_list[v]

            if len(neighbors) == 0:
                tau[v] = 0
                continue

            if k == 0:
                tau[v] = kCoreBaseStructuralDiversity._count_components_k0(neighbors, adj_sets)
            elif len(neighbors) < k:
                tau[v] = 0
            else:
                tau[v] = kCoreBaseStructuralDiversity._count_k_cores_optimized(neighbors, adj_sets, k)

        return tau

    @staticmethod
    def _count_components_k0(neighbors, adj_sets):
        if not neighbors:
            return 0

        neighbor_set = set(neighbors)
        visited = set()
        components = 0

        for u in neighbors:
            if u not in visited:

                components += 1
                queue = deque([u])
                visited.add(u)

                while queue:
                    curr = queue.popleft()

                    for next_node in adj_sets[curr]:
                        if next_node in neighbor_set and next_node not in visited:
                            visited.add(next_node)
                            queue.append(next_node)

        return components

    @staticmethod
    def _count_k_cores_fast(neighbors, adjecancy_sets, k):
        neighbor_set = set(neighbors)

        degrees = {}
        for u in neighbors:
            degrees[u] = len(adjecancy_sets[u] & neighbor_set)

        staying = set(neighbors)
        to_remove = deque([u for u in neighbors if degrees[u] < k])

        while to_remove:
            u = to_remove.popleft()
            if u not in staying:
                continue

            staying.remove(u)

            for v in adjecancy_sets[u]:
                if v in staying:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        to_remove.append(v)

        if not staying:
            return 0

        visited = set()
        components = 0

        for u in staying:
            if u not in visited:
                components += 1
                stack = [u]
                visited.add(u)

                while stack:
                    curr = stack.pop()
                    for next_node in adjecancy_sets[curr]:
                        if next_node in staying and next_node not in visited:
                            visited.add(next_node)
                            stack.append(next_node)

        return components

    @staticmethod
    def _count_k_cores_optimized(neighbors, adj_sets, k):
        if not neighbors:
            return 0

        neighbor_set = set(neighbors)
        n_neighbors = len(neighbors)

        neighbor_list = list(neighbors)
        neighbor_to_idx = {v: i for i, v in enumerate(neighbor_list)}

        degrees = [0] * n_neighbors
        for i, u in enumerate(neighbor_list):
            degrees[i] = len(adj_sets[u] & neighbor_set)

        active = [True] * n_neighbors
        queue = deque()

        for i in range(n_neighbors):
            if degrees[i] < k:
                queue.append(i)

        while queue:
            i = queue.popleft()
            if not active[i]:
                continue

            active[i] = False
            u = neighbor_list[i]

            for v in adj_sets[u]:
                if v in neighbor_set:
                    j = neighbor_to_idx[v]
                    if active[j]:
                        degrees[j] -= 1
                        if degrees[j] < k:
                            queue.append(j)

        visited = [False] * n_neighbors
        components = 0

        for i in range(n_neighbors):
            if active[i] and not visited[i]:
                components += 1
                stack = [i]
                visited[i] = True

                while stack:
                    curr_i = stack.pop()
                    curr_u = neighbor_list[curr_i]

                    for v in adj_sets[curr_u]:
                        if v in neighbor_set:
                            j = neighbor_to_idx[v]
                            if active[j] and not visited[j]:
                                visited[j] = True
                                stack.append(j)

        return components

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
