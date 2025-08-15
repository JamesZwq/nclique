#!/usr/bin/env python3
# Auto-generated for 5525875

STUDENT_ID = "5525875"
STUDENT_NAME = "Jingyu Sun"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque,defaultdict
################################################################################
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass
    @staticmethod
    def process(graph, k):
        # Prepare input arguments for each vertex
        tasks = [(node, graph.adj_list, k) for node in range(graph.vertex_num)]
        results = [kCoreBaseStructuralDiversity._evaluate(task) for task in tasks]

        return results

    @staticmethod
    def _evaluate(task):
        # Unpack task tuple
        node, adj_list, k = task
        neighbors = adj_list[node]

        # If the node has no neighbors, tau_k is zero
        if not neighbors:
            return 0

        # Create a subgraph induced by the node's neighbors
        neighbor_set = set(neighbors)
        subgraph = defaultdict(list)

        for u in neighbors:
            for v in adj_list[u]:
                if v in neighbor_set:
                    subgraph[u].append(v)

        # Count number of k-core components in the subgraph
        return kCoreBaseStructuralDiversity._count_k_core(subgraph, k)

    @staticmethod
    def _count_k_core(subgraph, k):
        # Return 0 if subgraph is empty
        if not subgraph:
            return 0

        # Initialize degree of each node
        degrees = {u: len(subgraph[u]) for u in subgraph}
        removed = set()
        queue = deque()

        # Mark nodes with degree less than k for removal
        for u in subgraph:
            if degrees[u] < k:
                removed.add(u)
                queue.append(u)

        # Remove nodes iteratively and update neighbors' degrees
        while queue:
            current = queue.popleft()
            for neighbor in subgraph[current]:
                if neighbor in removed:
                    continue
                degrees[neighbor] -= 1
                if degrees[neighbor] < k:
                    removed.add(neighbor)
                    queue.append(neighbor)

        # Nodes not removed form the remaining valid part of the subgraph
        remaining = [u for u in subgraph if u not in removed]
        if not remaining:
            return 0

        # Use BFS to count how many connected k-core components exist
        visited = set()
        count = 0

        for start in remaining:
            if start in visited:
                continue
            count += 1
            q = deque([start])
            visited.add(start)

            while q:
                u = q.popleft()
                for v in subgraph[u]:
                    if v not in removed and v not in visited:
                        visited.add(v)
                        q.append(v)

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
