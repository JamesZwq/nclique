#!/usr/bin/env python3
# Auto-generated for 5560088

STUDENT_ID = "5560088"
STUDENT_NAME = "Jiaao Yu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################
from collections import defaultdict, deque
class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        num_nodes = G.vertex_num
        result = [0] * num_nodes

        # go through every node
        for v in range(num_nodes):
            neighbors = G.adj_list[v]
            neighbor_set = set(neighbors)

            # build the edge list of neighbor-induced subgraph
            edge_list = []
            for u in neighbors:
                for w in G.adj_list[u]:
                  # avoid duplicate edges
                    if w in neighbor_set and u < w:
                        edge_list.append((u, w))

            # count k-core components in the neighbor subgraph
            result[v] = kCoreBaseStructuralDiversity.count_kcore_components(neighbor_set, edge_list, k)

        return result

    # count how many k-core connected components exist
    @staticmethod
    def count_kcore_components(nodes, edges, k):

        # adjacency list
        graph = defaultdict(set)
        for u, v in edges:
            graph[u].add(v)
            graph[v].add(u)

        degree = {node: len(graph[node]) for node in nodes}

        # remove < k node
        queue = deque([node for node in nodes if degree[node] < k])
        while queue:
            node = queue.popleft()
            for neighbor in graph[node]:
                graph[neighbor].remove(node)
                if degree[neighbor] >= k:
                    degree[neighbor] -= 1
                    if degree[neighbor] == k - 1:
                        queue.append(neighbor)
            del graph[node]
            del degree[node]

        # DFS
        visited = set()
        component_count = 0

        for node in graph:
            if node not in visited:
                component_count += 1
                stack = [node]
                visited.add(node)
                while stack:
                    current = stack.pop()
                    for neighbor in graph[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            stack.append(neighbor)

        return component_count
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
