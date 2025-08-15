#!/usr/bin/env python3
# Auto-generated for 5448611

STUDENT_ID = "5448611"
STUDENT_NAME = "Meizi Jiang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################


class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Calculates the k-core-based structural diversity using a highly efficient
        triangle-based pre-computation approach.

        Parameters
        ----------
        G : UndirectedUnweightedGraph
            The input graph.
        k : int
            The core number.

        Returns
        -------
        List[int]
            A list where the i-th element is the k-core-based structural
            diversity τ_k(i) for vertex i.
        """
        n = G.vertex_num
        adj = G.adj_list
        sd = [0] * n

        # --- Phase 1: Pre-computation of edge supports and triangles ---
        # support[(u, v)] stores the number of triangles containing edge (u, v).
        # triangles[(u, v)] stores the third vertices of these triangles.
        support = {}
        triangles = {}

        for u in range(n):
            for v in adj[u]:
                if u < v:
                    # Find common neighbors to count triangles for the edge (u, v)
                    # To optimize, iterate over the shorter neighbor list.
                    neighbors_u = set(adj[u])
                    neighbors_v = set(adj[v])
                    
                    if len(neighbors_u) < len(neighbors_v):
                        common_neighbors = [w for w in neighbors_u if w in neighbors_v]
                    else:
                        common_neighbors = [w for w in neighbors_v if w in neighbors_u]
                    
                    sup_count = len(common_neighbors)
                    support[(u, v)] = sup_count
                    support[(v, u)] = sup_count
                    
                    # Store the triangle info for efficient access later
                    triangles[(u, v)] = common_neighbors
                    triangles[(v, u)] = common_neighbors

        # --- Phase 2: Per-vertex calculation using pre-computed data ---
        for v in range(n):
            neighbors_of_v = adj[v]
            if not neighbors_of_v:
                sd[v] = 0
                continue

            # Special case for k=0: count connected components in nbr_v.
            # We can use the pre-computed triangles to define nbr_v's edges.
            if k == 0:
                visited = set()
                components = 0
                for node in neighbors_of_v:
                    if node not in visited:
                        components += 1
                        q = deque([node])
                        # BFS to visit connected component
                        visited.add(node)
                        while q:
                            curr = q.popleft()
                            # Neighbors of curr within nbr_v are defined by triangles with v
                            for neighbor in triangles.get((v, curr), []):
                                if neighbor in neighbors_of_v and neighbor not in visited:
                                    visited.add(neighbor)
                                    q.append(neighbor)
                sd[v] = components
                continue
            
            # General case for k > 0: Simulate k-core decomposition on nbr_v
            
            # 1. Get initial degrees from pre-computed support values.
            # This is a temporary copy for the local decomposition.
            temp_degrees = {u: support.get((v, u), 0) for u in neighbors_of_v}
            
            # 2. Initialize removal queue with nodes having degree < k.
            queue = deque([u for u, deg in temp_degrees.items() if deg < k])
            removed = set(queue)

            # 3. Iteratively remove nodes.
            while queue:
                u = queue.popleft()
                # When u is removed, the degrees of its neighbors in nbr_v decrease.
                # These neighbors are precisely the nodes w that form a triangle (v, u, w).
                for w in triangles.get((v, u), []):
                    if w in temp_degrees and w not in removed:
                        temp_degrees[w] -= 1
                        if temp_degrees[w] < k:
                            removed.add(w)
                            queue.append(w)
            
            # 4. Count connected components on the remaining nodes.
            remaining_nodes = set(neighbors_of_v) - removed
            if not remaining_nodes:
                sd[v] = 0
                continue

            visited = set()
            core_count = 0
            for node in remaining_nodes:
                if node not in visited:
                    core_count += 1
                    component_q = deque([node])
                    visited.add(node)
                    while component_q:
                        curr = component_q.popleft()
                        # Neighbors are defined by triangles with v, restricted to remaining nodes.
                        for neighbor in triangles.get((v, curr), []):
                            if neighbor in remaining_nodes and neighbor not in visited:
                                visited.add(neighbor)
                                component_q.append(neighbor)
            sd[v] = core_count

        return sd

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
