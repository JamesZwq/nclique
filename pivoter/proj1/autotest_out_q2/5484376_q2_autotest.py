#!/usr/bin/env python3
# Auto-generated for 5484376

STUDENT_ID = "5484376"
STUDENT_NAME = "Yuhang Liang"

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
            # Get neighbor-induced subgraph for vertex v
            neighbors = set(G.adj_list[v])
            if len(neighbors) == 0:
                sd[v] = 0
                continue

            # Create neighbor-induced subgraph
            neighbor_list = list(neighbors)
            neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbor_list)}

            # Build adjacency list for neighbor-induced subgraph
            sub_adj = [[] for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list):
                for neighbor in G.adj_list[u]:
                    if neighbor in neighbor_to_idx:
                        j = neighbor_to_idx[neighbor]
                        sub_adj[i].append(j)

            # Find all k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity._find_k_cores(sub_adj, k)
            sd[v] = len(k_cores)

        return sd

    @staticmethod
    def _find_k_cores(adj_list, k):
        """
        Find all maximal k-cores in the graph represented by adj_list
        Returns a list of k-cores, where each k-core is a set of vertices
        """
        n = len(adj_list)
        if n == 0:
            return []

        # Compute core numbers for all vertices
        core_numbers = kCoreBaseStructuralDiversity._compute_core_numbers(adj_list)

        # Find vertices with core number >= k
        k_core_vertices = set()
        for v in range(n):
            if core_numbers[v] >= k:
                k_core_vertices.add(v)

        if not k_core_vertices:
            return []

        # Find connected components in the k-core subgraph
        visited = [False] * n
        k_cores = []

        for v in k_core_vertices:
            if not visited[v]:
                # BFS to find connected component
                component = set()
                queue = deque([v])
                visited[v] = True

                while queue:
                    u = queue.popleft()
                    component.add(u)

                    for neighbor in adj_list[u]:
                        if neighbor in k_core_vertices and not visited[neighbor]:
                            visited[neighbor] = True
                            queue.append(neighbor)

                if component:
                    k_cores.append(component)

        return k_cores

    @staticmethod
    def _compute_core_numbers(adj_list):
        """
        Compute core numbers for all vertices using the standard algorithm
        """
        n = len(adj_list)
        degrees = [len(adj_list[v]) for v in range(n)]
        core_numbers = [0] * n

        # Create buckets for vertices by degree
        max_degree = max(degrees) if degrees else 0
        buckets = [[] for _ in range(max_degree + 1)]

        for v in range(n):
            buckets[degrees[v]].append(v)

        # Process vertices in order of increasing degree
        removed = [False] * n

        for d in range(max_degree + 1):
            while buckets[d]:
                v = buckets[d].pop()
                if removed[v]:
                    continue

                removed[v] = True
                core_numbers[v] = d

                # Update degrees of neighbors
                for neighbor in adj_list[v]:
                    if not removed[neighbor] and degrees[neighbor] > d:
                        # Remove from old bucket
                        if neighbor in buckets[degrees[neighbor]]:
                            buckets[degrees[neighbor]].remove(neighbor)

                        # Decrease degree and add to new bucket
                        degrees[neighbor] -= 1
                        if degrees[neighbor] >= d:
                            buckets[degrees[neighbor]].append(neighbor)

        return core_numbers


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
