#!/usr/bin/env python3
# Auto-generated for 5490499

STUDENT_ID = "5490499"
STUDENT_NAME = "David Rao"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def subgraph(G, neighbor):
        ################################################################################
        # Construct the induced subgraph for the neighbor set of a vertex.
        # Only keeps edges between neighbors.
        ################################################################################
        neighbor_set = set(neighbor)
        u_subgraph = defaultdict(set)
        for u in neighbor_set:
            for v in G.adj_list[u]:
                if v in neighbor_set:
                    u_subgraph[u].add(v)
        for u in neighbor_set:
            if u not in u_subgraph:
                u_subgraph[u] = []

        return u_subgraph
    
    @staticmethod
    def get_k_core(u_subgraph, k):
        ################################################################################
        # Extract the k-core of the subgraph.
        # Repeatedly remove nodes with degree less than k.
        ################################################################################
        degree = {}
        for u in u_subgraph:
            degree[u] = len(u_subgraph[u])
        
        queue = deque()
        for u in u_subgraph:
            if degree[u] < k:
                queue.append(u)
        
        set_subgraph = set(u_subgraph.keys())
        visited = set()

        while queue:
            u = queue.popleft()
            if u in visited:
                continue
            visited.add(u)
            for v in u_subgraph[u]:
                if v not in visited:
                    degree[v] -= 1
                    if degree[v] < k:
                        queue.append(v)
        set_subgraph = set_subgraph - visited
        
        if not set_subgraph:
            return set()
        
        return set_subgraph
    
    @staticmethod
    def count_vertex_tau(k_core_nodes, subgraph):
        ################################################################################
        # Count the number of connected components in the k-core subgraph using BFS.
        ################################################################################
        count = 0
        visited = set()
        for u in k_core_nodes:
            if u not in visited:
                count += 1
                visited.add(u)
                queue = deque([u])
                while queue:
                    v = queue.popleft()
                    for j in subgraph[v]:
                        if j not in visited and j in k_core_nodes:
                            visited.add(j)
                            queue.append(j)
        
        return count

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
            neighbor = list(G.adj_list[v])
            if not neighbor:
                sd[v] = 0
                continue

            # Build neighbor subgraph
            u_subgraph = kCoreBaseStructuralDiversity.subgraph(G, neighbor)

            # Find k-core nodes in subgraph
            k_core_nodes = kCoreBaseStructuralDiversity.get_k_core(u_subgraph, k)
            if not k_core_nodes:
                sd[v] = 0

            # Count components in the k-core
            count = kCoreBaseStructuralDiversity.count_vertex_tau(k_core_nodes, u_subgraph)

            sd[v] = count

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
