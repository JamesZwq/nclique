#!/usr/bin/env python3
# Auto-generated for 5486888

STUDENT_ID = "5486888"
STUDENT_NAME = "Lizhu Shen"

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
        if k < 0 or n == 0:
            return sd
        # compute core number for all vertice
        core_num = kCoreBaseStructuralDiversity._compute_core_num(G)
        # get vertices core number ≥ k+1
        k_plus_one_vertices = set(v for v in range(n) if core_num[v] >= k + 1)
        # if no (k+1)-core check k-core
        if not k_plus_one_vertices:
            return kCoreBaseStructuralDiversity._fallback_k_core_method(G, k, core_num)
        # find the connected component of (k+1)-core
        components = kCoreBaseStructuralDiversity._find_component(G, k_plus_one_vertices)
        # batch process each comp
        for component in components:
            kCoreBaseStructuralDiversity._process_comp_batch(G, component, k, sd)
        # fallback compute vertices not handled
        processed = set()
        for comp in components:
            for v in comp:
                for u in G.adj_list[v]:
                    processed.add(u)

        for v in range(n):
            if v not in processed or sd[v] == 0:
                sd[v] += kCoreBaseStructuralDiversity._compute_tau_k_local(G, v, k, core_num)

        return sd

    @staticmethod
    def _compute_core_num(G):
        """ optimize the cal of core numbers and use arrays instead of dictionaries to improve cache efficiency """
        n = G.vertex_num
        degree = [len(G.adj_list[v]) for v in range(n)]
        core_num = [0] * n
        # use count permutation optimize
        max_deg = max(degree, default=0)
        count = [0] * (max_deg + 1)
        for d in degree:
            count[d] += 1
        # set up
        start = [0] * (max_deg + 1)
        for i in range(1, max_deg + 1):
            start[i] = start[i - 1] + count[i - 1]
        # sort V
        sorted_ver = [0] * n
        pos = start[:]
        for v in range(n):
            sorted_ver[pos[degree[v]]] = v
            pos[degree[v]] += 1
        # peeling
        for i in range(n):
            v = sorted_ver[i]
            core_num[v] = degree[v]
            for u in G.adj_list[v]:
                if degree[u] > degree[v]:
                    degree[u] -= 1
        return core_num

    @staticmethod
    def _find_component(G, vertices):
        """ find com"""
        visited = set()
        components = []
        for v in vertices:
            if v not in visited:
                comp = []
                queue = deque([v])
                visited.add(v)
                while queue:
                    u = queue.popleft()
                    comp.append(u)
                    for w in G.adj_list[u]:
                        if w in vertices and w not in visited:
                            visited.add(w)
                            queue.append(w)
                components.append(comp)
        return components

    @staticmethod
    def _process_comp_batch(G, component, k, sd):
        """ batch process one (k+1) -core component """
        component_set = set(component)
        # collect all neighbors of the com
        neighbor_union = set()
        # for each
        for v in component:
            neighbor_union.update(G.adj_list[v])
        # compute how many k-core components its neighbor belong to
        for v in neighbor_union:
            nbrs = set(G.adj_list[v])
            component_neighbors = nbrs & component_set
            if len(component_neighbors) >= k:
                count = kCoreBaseStructuralDiversity._fast_k_core_count(G, component_neighbors, k, nbrs)
                sd[v] += count

    @staticmethod
    def _fast_k_core_count(G, kcore_candidates, k, permitted_neighbors):
        """ fast check before full peeling """
        if len(kcore_candidates) < k:
            return 0
        subgraph_degree = {}
        for v in kcore_candidates:
            deg = sum(1 for u in G.adj_list[v] if u in kcore_candidates and u in permitted_neighbors)
            if deg < k:
                return kCoreBaseStructuralDiversity._exact_k_core_count(G, kcore_candidates, k, permitted_neighbors)
            subgraph_degree[v] = deg
        # all degrees ≥ k skip peeling
        return kCoreBaseStructuralDiversity._count_components(G, kcore_candidates, permitted_neighbors)

    @staticmethod
    def _exact_k_core_count(G, vertices, k, permitted_neighbors):
        """ calculate the number of k-cores precisely """
        vertices = set(vertices)
        degree = {v: sum(1 for u in G.adj_list[v] if u in vertices and u in permitted_neighbors) for v in vertices}
        # peeling
        queue = deque(v for v in vertices if degree[v] < k)
        removed = set(queue)

        while queue:
            v = queue.popleft()
            for u in G.adj_list[v]:
                if u in vertices and u not in removed and u in permitted_neighbors:
                    degree[u] -= 1
                    if degree[u] < k:
                        queue.append(u)
                        removed.add(u)
        k_core_vertices = vertices - removed
        if not k_core_vertices:
            return 0
        return kCoreBaseStructuralDiversity._count_components(G, k_core_vertices, permitted_neighbors)

    @staticmethod
    def _count_components(G, vertices, permitted_neighbors):
        # count number of connected components
        visited = set()
        count = 0
        for v in vertices:
            if v not in visited:
                # bfs
                queue = deque([v])
                visited.add(v)
                while queue:
                    u = queue.popleft()
                    for w in G.adj_list[u]:
                        if w in vertices and w not in visited and w in permitted_neighbors:
                            queue.append(w)
                            visited.add(w)
                count += 1
        return count

    @staticmethod
    def _compute_tau_k_local(G, v, k, core_num):
        """ method for handling a single vertex """
        nbrs = [u for u in G.adj_list[v] if core_num[u] >= k]
        if len(nbrs) < k:
            return 0
        # create neighbors subgraph
        neighbor_set = set(nbrs)
        degree = {u: sum(1 for w in G.adj_list[u] if w in neighbor_set) for u in nbrs}
        # peeling
        queue = deque(u for u in nbrs if degree[u] < k)
        removed = set(queue)
        while queue:
            u = queue.popleft()
            for w in G.adj_list[u]:
                if w in neighbor_set and w not in removed:
                    degree[w] -= 1
                    if degree[w] < k:
                        queue.append(w)
                        removed.add(w)
        k_core_vertices = neighbor_set - removed
        if not k_core_vertices:
            return 0
        return kCoreBaseStructuralDiversity._count_components_fast(G, k_core_vertices)

    @staticmethod
    def _count_components_fast(G, vertices):
        """  count number of connected components in a set """
        if len(vertices) <= 1:
            return len(vertices)
        visited = set()
        count = 0
        # dfs
        for start in vertices:
            if start not in visited:
                # list instead of queue
                stack = [start]
                visited.add(start)
                while stack:
                    u = stack.pop()
                    for w in G.adj_list[u]:
                        if w in vertices and w not in visited:
                            stack.append(w)
                            visited.add(w)
                count += 1
        return count

    @staticmethod
    def _fallback_k_core_method(G, k, core_num):
        """ without (k+1)-core fallback method """
        n = G.vertex_num
        sd = [0] * n
        for v in range(n):
            nbrs = [u for u in G.adj_list[v] if core_num[u] >= k]
            if len(nbrs) >= k:
                sd[v] = kCoreBaseStructuralDiversity._compute_tau_k_local(G, v, k, core_num)
        return sd

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
