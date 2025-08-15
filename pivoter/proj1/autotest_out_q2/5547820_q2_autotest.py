#!/usr/bin/env python3
# Auto-generated for 5547820

STUDENT_ID = "5547820"
STUDENT_NAME = "Xuanling Wang"

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

        # Get neighbors of all vertex, and construct them to subgraphs
        for v in range(n):
            v_neighbors = G.adj_list[v]
            # If there is no neighbor
            if not v_neighbors:
                sd[v] = 0
                continue
            # vertice v is not included in subgraph, we just construct a adjlist here
            subGraph = {u: set() for u in v_neighbors}

            #
            for u in v_neighbors:
                for u_neighbors in G.adj_list[u]:
                    if u_neighbors in v_neighbors:
                        subGraph[u].add(u_neighbors)
            
            # compute K-core if there is valid K-core and then decomposite to count the number of K-core
            sd[v] = kCoreBaseStructuralDiversity.computeAndCountKcore(subGraph, k)

        return sd

    @staticmethod
    def computeAndCountKcore(subGraph, k):
        # Store the vertices that will be removed during k-core peeling
        deletedVertex = set()
        queue = deque()
        # Compute the initial degree of each vertice in the subgraph, and Extract K-core vertex
        degrees = {u: len(neighbors) for u, neighbors in subGraph.items()}
        
        # Initialize the peeling process by adding all nodes with degree < k into the queue
        for u in subGraph:
            if degrees[u] < k:
                deletedVertex.add(u)
                queue.append(u)
        
        # Iteratively remove vertex with degree less than k 
        # When a vertice is removed, its neighbors' degrees are updated
        # If the degree of a neighbor drops below k, it is also removed
        while queue:
            u = queue.popleft()
            for v in subGraph[u]:
                if v not in deletedVertex:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        deletedVertex.add(v)
                        queue.append(v)
        
        # After peeling, the remained vertex are those in the k-core
        # These will be used to count the number of connected components
        remainedVertex = set(subGraph.keys()) - deletedVertex

        return kCoreBaseStructuralDiversity.countKcore(subGraph, remainedVertex)
    
    @staticmethod
    def countKcore(subGraph, remainedVertex):
        # Track visited nodes to avoid counting the same component multiple times
        visited = set()
        count = 0

        # Perform BFS from each unvisited node to count connected components
        # Each BFS traversal corresponds to one distinct k-core component
        for u in remainedVertex:
            if u not in visited:
                count += 1
                # Initialize BFS from this node
                queue = deque()
                queue.append(u)
                visited.add(u)

                # Standard BFS to traverse the entire connected component
                # Only vertex in remainedVertex are considered valid
                while queue:
                    v = queue.popleft()
                    for w in subGraph[v]:
                        if w in remainedVertex and w not in visited:
                            visited.add(w)
                            queue.append(w)
        
        return count

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
