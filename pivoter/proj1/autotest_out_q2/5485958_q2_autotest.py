#!/usr/bin/env python3
# Auto-generated for 5485958

STUDENT_ID = "5485958"
STUDENT_NAME = "Haoyan Wei"

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
        
        if k == 0:
            return [0] * n
        
        # Preprocess: mark vertices with global degree < k and Calculate global degree for each vertex
        global_deg = [len(adj) for adj in G.adj_list]
        low_degree_nodes = [deg < k for deg in global_deg]
        
        for v in range(n):
            # If the degree of vertex v is less than k, then τ_k(v) = 0
            if global_deg[v] < k:
                sd[v] = 0
                continue
            
            # Get candidates with global degree >= k
            candidates = [u for u in G.adj_list[v] if not low_degree_nodes[u]]
            if not candidates:
                sd[v] = 0
                continue
                
            candidate_set = set(candidates)
            deg_in_subgraph = {}
            queue = deque()
            
            # Initialize the degree of each candidate in the subgraph
            for u in candidates:
                # Calculate the degree of u in the induced subgraph
                cnt = 0
                for w in G.adj_list[u]:
                    if w in candidate_set:
                        cnt += 1
                deg_in_subgraph[u] = cnt
                if cnt < k:
                    queue.append(u)
            
            # K-core decomposition
            deleted = set()
            while queue:
                u = queue.popleft()
                if u in deleted:
                    continue
                deleted.add(u)
                for w in G.adj_list[u]:
                    if w in candidate_set and w not in deleted:
                        deg_in_subgraph[w] -= 1
                        if deg_in_subgraph[w] < k:
                            queue.append(w)
            
            # Get the remaining nodes after deletion
            remaining = candidate_set - deleted
            if not remaining:
                sd[v] = 0
                continue
                
            # Calculate the number of connected components in the remaining subgraph
            visited = {node: False for node in remaining}
            comp_count = 0
            for node in remaining:
                if not visited[node]:
                    comp_count += 1
                    stack = [node]
                    visited[node] = True
                    while stack:
                        cur = stack.pop()
                        for neighbor in G.adj_list[cur]:
                            if neighbor in remaining and not visited[neighbor]:
                                visited[neighbor] = True
                                stack.append(neighbor)
            sd[v] = comp_count
        
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
