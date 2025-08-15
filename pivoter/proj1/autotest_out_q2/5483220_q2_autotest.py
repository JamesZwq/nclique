#!/usr/bin/env python3
# Auto-generated for 5483220

STUDENT_ID = "5483220"
STUDENT_NAME = "Xin Xie"

# ======= 学生代码 =======
from collections import deque

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        # each vertex resut
        result = []

        # iter each vertex
        for i in range(G.vertex_num):
            # get all neigh of i
            my_neighbors = G.adj_list[i]

            # if no neigh，the number of k-core is 0
            if len(my_neighbors) == 0:
                result.append(0)
                continue

            # build new number of neigh vertex
            old_to_new = {}  # old -> new
            new_to_old = {}  # new -> old
            for idx, neighbor in enumerate(my_neighbors):
                old_to_new[neighbor] = idx
                new_to_old[idx] = neighbor

            # build graph
            neighbor_graph = [[] for _ in range(len(my_neighbors))]

            # check if all neigh has edge
            for neighbor in my_neighbors:
                new_id = old_to_new[neighbor]
                # check all neigh node of this node
                for adj_vertex in G.adj_list[neighbor]:
                    # if the adj node is in old to new
                    if adj_vertex in old_to_new:
                        adj_new_id = old_to_new[adj_vertex]
                        neighbor_graph[new_id].append(adj_new_id)

            # count degree of neigh subgraph vertext
            vertex_degrees = [0] * len(neighbor_graph)
            for j in range(len(neighbor_graph)):
                vertex_degrees[j] = len(neighbor_graph[j])

            # store removed node
            removed_vertices = [False] * len(neighbor_graph)

            # if degree < k, then remove
            changed = True
            while changed:
                changed = False
                for j in range(len(neighbor_graph)):
                    if not removed_vertices[j] and vertex_degrees[j] < k:
                        # remove the vertex
                        removed_vertices[j] = True
                        changed = True
                        # reduce its neigh degree
                        for neighbor_idx in neighbor_graph[j]:
                            vertex_degrees[neighbor_idx] -= 1

            # count remaining connected componet
            visited_vertices = [False] * len(neighbor_graph)
            connected_components = 0

            for j in range(len(neighbor_graph)):
                if not removed_vertices[j] and not visited_vertices[j]:
                    # fing a new, use BFS
                    connected_components += 1
                    bfs_queue = deque()
                    bfs_queue.append(j)
                    visited_vertices[j] = True

                    while bfs_queue:
                        current_vertex = bfs_queue.popleft()
                        for neighbor_vertex in neighbor_graph[current_vertex]:
                            if not removed_vertices[neighbor_vertex] and not visited_vertices[neighbor_vertex]:
                                visited_vertices[neighbor_vertex] = True
                                bfs_queue.append(neighbor_vertex)

            result.append(connected_components)

        return result

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
