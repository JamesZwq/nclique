#!/usr/bin/env python3
# Auto-generated for 5547616

STUDENT_ID = "5547616"
STUDENT_NAME = "Mingyuan Sun"

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
            # Get the neighbors of vertex v
            neighbors = kCoreBaseStructuralDiversity._get_neighbors(G, v)
            if not neighbors:
                continue
            # Build the sub - graph based on the neighbors of vertex v
            sub_graph = kCoreBaseStructuralDiversity._build_sub_graph(G, neighbors)
            # Calculate the number of k - cores in the sub - graph
            sd[v] = kCoreBaseStructuralDiversity._count_k_core(sub_graph, k)

        return sd

    @staticmethod
    def _get_neighbors(G, v):
        #Get the set of neighbors of vertex v
        neighbors = G.adj_list[v]
        if not neighbors:
            return set()
        return set(neighbors)

    @staticmethod
    def _build_sub_graph(G, neighbors):
        #Build a sub - graph based on the given set of neighbors
        sub_G = {u: [] for u in neighbors}
        for u in neighbors:
            for neighbor_of_u in G.adj_list[u]:
                if neighbor_of_u in neighbors:
                    sub_G[u].append(neighbor_of_u)
        return sub_G

    @staticmethod
    def _count_k_core(sub_G, k):
        #Calculate the number of k - cores in the sub - graph
        if not sub_G:
            return 0
        # Compute the degree of each vertex in the sub - graph
        degrees = kCoreBaseStructuralDiversity._compute_degrees(sub_G)
        # Find the vertices whose degree is less than k
        vertices_to_delete = kCoreBaseStructuralDiversity._find_vertices_to_delete(degrees, k)
        # Iteratively delete vertices whose degree is less than k
        deleted = kCoreBaseStructuralDiversity._delete_vertices(sub_G, degrees, vertices_to_delete, k)
        # Get the remaining vertices after deletion
        remains = kCoreBaseStructuralDiversity._get_remaining_vertices(sub_G, deleted)
        if not remains:
            return 0
        # Count the number of connected components of the remaining vertices
        return kCoreBaseStructuralDiversity._count_connected_components(sub_G, remains, deleted)

    @staticmethod
    def _compute_degrees(sub_G):
        #Compute the degree of each vertex in the subgraph
        return {u: len(neighbors) for u, neighbors in sub_G.items()}

    @staticmethod
    def _find_vertices_to_delete(degrees, k):
        #Find the vertices whose degree is less than k
        vertices = deque()
        deleted = set()
        for u, degree in degrees.items():
            if degree < k:
                vertices.append(u)
                deleted.add(u)
        return vertices, deleted

    @staticmethod
    def _delete_vertices(sub_G, degrees, vertices, k):
        #Iteratively delete vertices whose degree is less than k
        deleted = vertices[1]
        vertices = vertices[0]
        while vertices:
            u = vertices.popleft()
            for v in sub_G[u]:
                if v not in deleted:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        deleted.add(v)
                        vertices.append(v)
        return deleted

    @staticmethod
    def _get_remaining_vertices(sub_G, deleted):
        #Get the remaining vertices after deletion
        return [u for u in sub_G if u not in deleted]

    @staticmethod
    def _count_connected_components(sub_G, remains, deleted):
        #Count the number of connected components of the remaining vertices
        visited = set()
        k_core_count = 0
        for node in remains:
            if node not in visited:
                k_core_count += 1
                queue = deque([node])
                visited.add(node)
                while queue:
                    u = queue.popleft()
                    for v in sub_G[u]:
                        if v not in deleted and v not in visited:
                            visited.add(v)
                            queue.append(v)
        return k_core_count

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
