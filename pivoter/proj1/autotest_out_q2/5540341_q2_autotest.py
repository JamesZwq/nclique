#!/usr/bin/env python3
# Auto-generated for 5540341

STUDENT_ID = "5540341"
STUDENT_NAME = "Runfeng Huang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
from typing import List
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k: int) -> List[int]:
        # Calculate the global core number of each vertex in the graph, that is, which layer of k-core the node is in.
        def _global_core_numbers(graph):
            n = graph.vertex_num
            deg_list = [len(graph.adj_list[i]) for i in range(n)]
            core_val = [0] * n
            max_deg = max(deg_list) if n > 0 else 0

            freq = [0] * (max_deg + 1)
            for d in deg_list:
                freq[d] += 1

            start_idx = [0] * (max_deg + 1)
            for i in range(1, max_deg + 1):
                start_idx[i] = start_idx[i - 1] + freq[i - 1]

            vert_seq = [0] * n
            pos_in_seq = [0] * n
            for v in range(n):
                d = deg_list[v]
                vert_seq[start_idx[d]] = v
                pos_in_seq[v] = start_idx[d]
                start_idx[d] += 1
            for i in range(max_deg, 0, -1):
                start_idx[i] = start_idx[i - 1]
            start_idx[0] = 0

            for i in range(n):
                v = vert_seq[i]
                for u in graph.adj_list[v]:
                    if deg_list[u] > deg_list[v]:
                        du = deg_list[u]
                        pu, pw = pos_in_seq[u], start_idx[du]
                        w = vert_seq[pw]
                        if u != w:
                            pos_in_seq[u], pos_in_seq[w] = pw, pu
                            vert_seq[pu], vert_seq[pw] = w, u
                        start_idx[du] += 1
                        deg_list[u] -= 1
                core_val[v] = deg_list[v]
            return core_val

        # Given a set of vertices, find all connected components.
        def _bfs_components(graph, node_group):
            visited = set()
            comps = []
            for v in node_group:
                if v not in visited:
                    part = []
                    dq = deque([v])
                    visited.add(v)
                    while dq:
                        u = dq.popleft()
                        part.append(u)
                        for nb in graph.adj_list[u]:
                            if nb in node_group and nb not in visited:
                                visited.add(nb)
                                dq.append(nb)
                    comps.append(part)
            return comps
        # Count the number of connected components in a subgraph (consisting of nodeset).
        def _dfs_component_count(graph, nodeset):
            visited = set()
            count = 0
            for v in nodeset:
                if v not in visited:
                    stack = [v]
                    visited.add(v)
                    while stack:
                        u = stack.pop()
                        for w in graph.adj_list[u]:
                            if w in nodeset and w not in visited:
                                visited.add(w)
                                stack.append(w)
                    count += 1
            return count

        # In a neighbor set, count the number of connected components in its k-core.
        def _local_k_core_count(graph, subset, kval):
            if len(subset) < kval:
                return 0
            s = set(subset)
            degree = {v: sum(1 for u in graph.adj_list[v] if u in s) for v in s}
            removed = set()
            dq = deque([v for v in s if degree[v] < kval])
            removed.update(dq)

            while dq:
                v = dq.popleft()
                for u in graph.adj_list[v]:
                    if u in s and u not in removed:
                        degree[u] -= 1
                        if degree[u] < kval:
                            dq.append(u)
                            removed.add(u)
            remain = s - removed
            return _dfs_component_count(graph, remain)

        # For all connected components of the (k+1)-core, the structural diversity of each relevant node is calculated in batches.
        def _process_comps(graph, component_list, core_array, kval, diversity):
            for part in component_list:
                related = set(part)
                for node in part:
                    related.update(graph.adj_list[node])
                for v in related:
                    neighbors = [u for u in graph.adj_list[v] if u in related and core_array[u] >= kval]
                    if len(neighbors) < kval:
                        continue
                    c = _local_k_core_count(graph, neighbors, kval)
                    diversity[v] = max(diversity[v], c)
        # Basic validity check (empty graph or invalid k)
        num_nodes = G.vertex_num
        if num_nodes == 0 or k < 0:
            return [0] * num_nodes
        # Initialize the return array + global core value calculation
        tau_k = [0] * num_nodes
        core_arr = _global_core_numbers(G)
        # Find all nodes belonging to the (k+1)-core
        # If (k+1)-core is empty, all diversity is 0
        S_k1 = [v for v in range(num_nodes) if core_arr[v] >= k + 1]
        if not S_k1:
            return tau_k

        # Each (k+1)-core connected component is batch processed and the structural diversity value is finally returned.
        comps_k1 = _bfs_components(G, S_k1)
        _process_comps(G, comps_k1, core_arr, k, tau_k)
        return tau_k


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
