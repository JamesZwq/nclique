#!/usr/bin/env python3
# Auto-generated for 5505772

STUDENT_ID = "5505772"
STUDENT_NAME = "Jianzhang Wu"

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
        n = G.vertex_num
        result = [0] * n
        for vertex in range(n):
            neighbor_set = set(G.adj_list[vertex])        
            # Construct the neighbor-induced subgraph and mark the valid_nodes
            subgraph = kCoreBaseStructuralDiversity.buildSubgraph(G, neighbor_set)
            # Find all k-core connected components in the subgraph
            k_cores = kCoreBaseStructuralDiversity.kCores(subgraph, k)
            result[vertex] = len(k_cores)

        return result

    @staticmethod
    def buildSubgraph(G, nodes):
        edges = []
        # Collect all undirected edges whose endpoints lie within the nodes
        for node in nodes:
            for neighbor  in G.adj_list[node]:
                if neighbor in nodes and node < neighbor:
                    edges.append((node, neighbor))
        m = len(edges)
        edge_list = [(G.vertex_num, m)] + edges
        # Construct the subgraph and mark its valid_nodes
        subgraph = UndirectedUnweightedGraph(edge_list)
        subgraph.valid_nodes = set(nodes)
        return subgraph

    # Extract all remaining connected components in the subgraph G, and return them as a list
    @staticmethod
    def kCores(G, k):
        if hasattr(G, 'valid_nodes'):
            nodes = G.valid_nodes
        else:
            nodes = set(range(G.vertex_num))

        # Initialize the degree of each node
        degree = {}
        for vertex in nodes:
            count = 0
            for neighbor in G.adj_list[vertex]:
                if neighbor in nodes:
                    count += 1
            degree[vertex] = count
            
        removed = {}
        for vertex in nodes:
            removed[vertex] = False

        # Initialize a queue and enqueue all nodes whose degree is less than k
        queue = deque()
        for vertex in nodes:
            if degree[vertex] < k:
                queue.append(vertex)


        # Iteratively remove nodes with degree < k and update the degrees of their neighbors
        while len(queue) > 0:
            current = queue.popleft()
            if removed[current]:
                continue
            removed[current] = True
            # Update the degree of each neighbor
            for neighbor in G.adj_list[current]:
                if neighbor in nodes and not removed[neighbor]:
                    degree[neighbor] -= 1
                    if degree[neighbor] < k:
                        queue.append(neighbor)
        
        visited = {}
        for vertex in nodes:
            visited[vertex] = False
        components = []
    
        def bfs(start):
            component = set()
            queue = deque()
            queue.append(start)
            visited[start] = True
        
            while len(queue) > 0:
                current = queue.popleft()
                component.add(current)
                for neighbor in G.adj_list[current]:
                    if neighbor in nodes and not removed[neighbor] and not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)
            return component
            
        # For each node that is neither removed nor visited, initiate a BFS to extract its connected component
        for vertex in nodes:
            if not removed[vertex] and not visited[vertex]:
                component = bfs(vertex)
                components.append(component)

        return components


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
