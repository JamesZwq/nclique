#!/usr/bin/env python3
# Auto-generated for 5461235

STUDENT_ID = "5461235"
STUDENT_NAME = "Haozhou Ni"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################
from collections import deque, defaultdict
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
            neighbors = G.adj_list[v]
            if not neighbors:
                continue

            # Step 1: Build Neighbor Induced Subgraph
            nbr_nodes = set(neighbors)
            nbr_adj = defaultdict(list)
            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in nbr_nodes:
                        nbr_adj[u].append(w)

            # Step 2: Get the connected component of the Induced Subgraph
            connected_components = kCoreBaseStructuralDiversity.get_connected_components(nbr_nodes, nbr_adj)

            # Step 3: Get k-core for each connected component
            count = 0
            for i in connected_components:
                # Build subgraph adjacency list
                i_adj = {}
                for u in i:
                    i_adj[u] = []
                    for w in nbr_adj[u]:
                        if w in i:
                            i_adj[u].append(w)
                
                core = kCoreBaseStructuralDiversity.get_k_core(i_adj, k)
                if core:
                    # Repartition the connected components within the k-core
                    core_adj = {}
                    for u in core:
                        core_adj[u] = []
                        for w in i_adj[u]:
                            if w in core:
                                core_adj[u].append(w)
                    connected_cores = kCoreBaseStructuralDiversity.get_connected_components(core, core_adj)
                    count += len(connected_cores)

            sd[v] = count

        return sd
    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def get_connected_components(nodes, adj):
        """BFS Find the connected component"""
        visited = set()
        components = []

        for node in nodes:
            if node not in visited:
                comp = set()
                queue = deque([node])
                visited.add(node)
                while queue:
                    u = queue.popleft()
                    comp.add(u)
                    for v in adj[u]:
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)
                components.append(comp)
        return components
    
    @staticmethod
    def get_k_core(adj, k):
        """Extract k-core from subgraph"""
        degree = {u: len(adj[u]) for u in adj}
        q = deque([u for u in adj if degree[u] < k])
        removed = set()

        while q:
            u = q.popleft()
            removed.add(u)
            for v in adj[u]:
                if v not in removed:
                    degree[v] -= 1
                    if degree[v] == k - 1:
                        q.append(v)
        
        return set(adj.keys()) - removed

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
