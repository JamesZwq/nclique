#!/usr/bin/env python3
# Auto-generated for 5510492

STUDENT_ID = "5510492"
STUDENT_NAME = "Kaiting Yang"

# ======= 学生代码 =======
from collections import deque
from typing import List, Dict

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k: int) -> List[int]:
        n = G.vertex_num
        tau = [0] * n  

        for v in range(n):
            neighbours = G.adj_list[v]  

            if not neighbours:
                tau[v] = 0
                continue

            # Construct the neighbor-induced subgraph 
            # This subgraph includes only the nodes in the neighbor set and the connections between them
            sub_G = {u: [] for u in neighbours}  
            for u in neighbours:
                for w in G.adj_list[u]:  
                    if w in neighbours:  
                        sub_G[u].append(w)

            # Call the function to compute the number of k-core connected components in the subgraph
            tau[v] = kCoreBaseStructuralDiversity._compute_k_core(sub_G, k)
        return tau

    @staticmethod
    def _compute_k_core(sub_G: Dict[int, List[int]], k: int) -> int:
        """
        Subfunction: Compute the number of connected components in the given subgraph that belong to the k-core
        Parameters:
            sub_G: The adjacency list of the subgraph (constructed only from the neighbors of a single node)
            k: Threshold for the k-core
        Returns:
            The number of connected components in sub_G that are part of the k-core (τ_k(v))
        """

        if not sub_G:
            return 0

        # Initialize the degree of each node
        degrees = {u: len(neigh) for u, neigh in sub_G.items()}
        deleted = set() 
        queue = deque()

        # add all nodes with degree < k to the removal queue
        for u in sub_G:
            if degrees[u] < k:
                queue.append(u)
                deleted.add(u)

        # Use BFS to continuously remove nodes with degree < k.
        # Each time a node is removed, update the degrees of its neighbors and check if their degree also falls below k.
        while queue:
            u = queue.popleft()
            for v in sub_G.get(u, []):  
                if v not in deleted:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        deleted.add(v)
                        queue.append(v)

        # Only the remaining nodes are part of the k-core
        remains = [u for u in sub_G if u not in deleted]
        if not remains:
            return 0  

        # Perform BFS traversal on the nodes in the k-core to count the number of connected components
        visited = set()
        count = 0

        for u in remains:
            if u not in visited:
                count += 1  
                queue = deque([u])
                visited.add(u)

                while queue:
                    curr = queue.popleft()
                    for nei in sub_G.get(curr, []):
                        if nei in remains and nei not in visited:
                            visited.add(nei)
                            queue.append(nei)
        return count

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
