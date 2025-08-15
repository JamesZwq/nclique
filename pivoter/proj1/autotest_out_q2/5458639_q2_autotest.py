#!/usr/bin/env python3
# Auto-generated for 5458639

STUDENT_ID = "5458639"
STUDENT_NAME = "Yuxuan Du"

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
        result = [0] * n

        # For each vertex v, compute τ_k(v) = |K_k(nbr_v)|
        # where nbr_v is the neighbor-induced subgraph of v
        # and K_k(nbr_v) is the set of k-cores in nbr_v

        for v in range(n):
            neighbors = G.adj_list[v]
            if len(neighbors) == 0:
                result[v] = 0
                continue

            # Build neighbor-induced subgraph
            neighbor_set = set(neighbors)
            neighbor_list = list(neighbors)
            neighbor_to_idx = {neighbor: i for i, neighbor in enumerate(neighbor_list)}

            # Create adjacency list for neighbor-induced subgraph
            nbr_adj = [[] for _ in range(len(neighbor_list))]
            for i, u in enumerate(neighbor_list):
                for neighbor_of_u in G.adj_list[u]:
                    if neighbor_of_u in neighbor_set:
                        j = neighbor_to_idx[neighbor_of_u]
                        nbr_adj[i].append(j)

            # Find k-cores in the neighbor-induced subgraph
            k_cores = kCoreBaseStructuralDiversity._find_all_k_cores(nbr_adj, k)
            result[v] = len(k_cores)

        return result

    @staticmethod
    def _find_all_k_cores(adj_list, k):
        """
        Find all k-cores in a graph represented by adjacency list
        Returns a list of k-core components (each component is a set of vertices)
        """
        n = len(adj_list)
        if n == 0:
            return []

        # Step 1: Find vertices that belong to some k-core
        k_core_vertices = kCoreBaseStructuralDiversity._find_k_core_vertices(adj_list, k)

        if not k_core_vertices:
            return []

        # Step 2: Find connected components within k-core vertices
        visited = [False] * n
        components = []

        for v in k_core_vertices:
            if not visited[v]:
                # BFS to find connected component
                component = set()
                queue = deque([v])
                visited[v] = True

                while queue:
                    curr = queue.popleft()
                    component.add(curr)

                    for neighbor in adj_list[curr]:
                        if not visited[neighbor] and neighbor in k_core_vertices:
                            visited[neighbor] = True
                            queue.append(neighbor)

                if component:
                    components.append(component)

        return components

    @staticmethod
    def _find_k_core_vertices(adj_list, k):
        """
        Find all vertices that belong to k-core using iterative removal
        """
        n = len(adj_list)
        degrees = [len(adj_list[v]) for v in range(n)]
        removed = [False] * n
        queue = deque()

        # Initially add all vertices with degree < k
        for v in range(n):
            if degrees[v] < k:
                queue.append(v)
                removed[v] = True

        # Iteratively remove vertices
        while queue:
            v = queue.popleft()

            # Update degrees of neighbors
            for neighbor in adj_list[v]:
                if not removed[neighbor]:
                    degrees[neighbor] -= 1
                    if degrees[neighbor] < k:
                        queue.append(neighbor)
                        removed[neighbor] = True

        # Collect remaining vertices (k-core vertices)
        return set(v for v in range(n) if not removed[v])

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
