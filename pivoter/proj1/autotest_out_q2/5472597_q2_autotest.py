#!/usr/bin/env python3
# Auto-generated for 5472597

STUDENT_ID = "5472597"
STUDENT_NAME = "Yue Liu"

# ======= 学生代码 =======
import collections

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def process(G, k):
        """
        Compute structural diversity τₖ(v) for each node based on K-core

        """
        n = G.vertex_num
        sd = [0] * n  # Store structural diversity results

        # Iterate all nodes
        for v in range(n):
            # Get neighbors of node v
            neighbors = G.adj_list[v]

            # If no neighbors, τₖ(v) = 0
            if not neighbors:
                sd[v] = 0
                continue
            # Build neighbor-induced subgraph
            # Use original ID to subgraph index mapping
            vertex_map = {u: idx for idx, u in enumerate(neighbors)}

            #  Length is #neighbors; index is subgraph index
            subgraph_adj = [[] for _ in range(len(neighbors))]

            # Fill subgraph adjacency list
            # Build edges among neighbors
            for idx_u, u in enumerate(neighbors): # Iterate neighbors; idx_u is index of u
                # Only scan u's neighbors in original graph
                for w in G.adj_list[u]:
                    # Check if w is also v’s neighbor
                    # Use u < w to avoid duplicates (undirected graph)
                    if w in vertex_map and u < w:
                        idx_w = vertex_map[w]
                        subgraph_adj[idx_u].append(idx_w)
                        subgraph_adj[idx_w].append(idx_u)


            # Run K-core on neighbor-induced subgraph
            degrees = [len(adj) for adj in subgraph_adj]
            removed = [False] * len(degrees)
            q = collections.deque()

            for i, deg in enumerate(degrees):
                if deg < k:
                    q.append(i)
                    removed[i] = True

            while q:
                node_idx = q.popleft() #  Run K-core on neighbor-induced subgraph
                for neighbor_idx in subgraph_adj[node_idx]: # removed[] indexed same as subgraph_adj
                    if not removed[neighbor_idx]:
                        degrees[neighbor_idx] -= 1
                        if degrees[neighbor_idx] < k:
                            removed[neighbor_idx] = True
                            q.append(neighbor_idx)

            # Count connected components in K-core
            # Map old subgraph index to new kcore index
            kcore_internal_map = {} # old subgraph index → new k-core index
            current_new_idx = 0
            for i in range(len(degrees)): # Iterate all original subgraph indices
                if not removed[i]:
                    kcore_internal_map[i] = current_new_idx
                    current_new_idx += 1


            if not kcore_internal_map: # If no remaining nodes (empty map)
                sd[v] = 0
                continue

            # Build K-core adjacency with new indices
            kcore_adj = [[] for _ in range(len(kcore_internal_map))] # Length = k-core nodes

            for original_subgraph_idx in kcore_internal_map.keys(): # Loop over original indices in K-core
                new_idx_u = kcore_internal_map[original_subgraph_idx] # Get new index in kcore_adj
                for neighbor_original_subgraph_idx in subgraph_adj[original_subgraph_idx]: # Iterate its neighbors in original subgraph
                    # Check if neighbor is still in K-core
                    if not removed[neighbor_original_subgraph_idx]:
                        new_idx_w = kcore_internal_map[neighbor_original_subgraph_idx] # Get new index of neighbor
                        kcore_adj[new_idx_u].append(new_idx_w)

            # Count components using new indices
            visited = [False] * len(kcore_internal_map) # visited[] size = #k-core nodes
            component_count = 0

            # Iterate over new k-core indices
            for i in range(len(kcore_internal_map)):
                if not visited[i]:
                    component_count += 1
                    queue = collections.deque([i]) # Push new index to queue
                    visited[i] = True

                    while queue:
                        cur_new_idx = queue.popleft()
                        for neighbor_new_idx in kcore_adj[cur_new_idx]:
                            # neighbor_new_idx is already correct k-core index
                            if not visited[neighbor_new_idx]:
                                visited[neighbor_new_idx] = True
                                queue.append(neighbor_new_idx)

            sd[v] = component_count

        return sd

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
