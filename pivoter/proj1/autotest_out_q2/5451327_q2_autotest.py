#!/usr/bin/env python3
# Auto-generated for 5451327

STUDENT_ID = "5451327"
STUDENT_NAME = "Zhuoyi Xiao"

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
        diversity = [0] * n

        # Precompute neighbor sets so we don't keep rebuilding them
        nbr_sets = [None] * n
        for v in range(n):
            nbr_sets[v] = set(G.adj_list[v])

        for v in range(n):
            neigh_set = nbr_sets[v]
            if not neigh_set:
                # No neighbors
                continue

            # Induced subgraph on v's neighbors
            subgraph = kCoreBaseStructuralDiversity.build_induced_subgraph(nbr_sets, neigh_set)

            # Quick exit for k = 0
            if k == 0:
                diversity[v] = kCoreBaseStructuralDiversity.count_components(subgraph)
                continue

            # Peel to k-core
            kept_nodes = kCoreBaseStructuralDiversity.peel_to_k_core(subgraph, k)
            if not kept_nodes:
                # Everything got peeled off
                continue

            # Filter subgraph to nodes that remain, then count CCs
            filtered = {}
            for u, neighs in subgraph.items():
                if u in kept_nodes:
                    kept_list = []
                    for w in neighs:
                        if w in kept_nodes:
                            kept_list.append(w)
                    filtered[u] = kept_list
            diversity[v] = kCoreBaseStructuralDiversity.count_components(filtered)
        return diversity


    @staticmethod
    def build_induced_subgraph(nbr_sets, nodes):
        """
        Build the induced subgraph over 'nodes' using the global neighbor sets.
        Returns: dict(node -> list_of_neighbors_within_nodes)
        """
        induced = {}
        for u in nodes:
            filtered_list = []
            for w in nbr_sets[u]:
                if w in nodes:
                    filtered_list.append(w)
            induced[u] = filtered_list
        return induced

    @staticmethod
    def peel_to_k_core(adj, k):
        """
        Standard k-core peeling:
          - start with all nodes whose degree < k
          - iteratively remove them and update others' degrees
        Returns the set of nodes that remain after peeling.
        """
        deg_map = kCoreBaseStructuralDiversity.build_degree_map(adj)
        queue, removed = kCoreBaseStructuralDiversity.init_prune_queue(deg_map, k)
        kCoreBaseStructuralDiversity.prune_nodes(adj, deg_map, queue, removed, k)
        return kCoreBaseStructuralDiversity.collect_remaining_nodes(adj, removed)

    @staticmethod
    def build_degree_map(adj):
        """Create a simple node->degree dict."""
        deg_map = {}
        for u, neighs in adj.items():
            deg_map[u] = len(neighs)
        return deg_map

    @staticmethod
    def init_prune_queue(deg_map, k):
        """
        Seed the queue with nodes whose degree is already < k.
        Also mark them as removed.
        """
        queue = deque()
        removed = set()
        for u, d in deg_map.items():
            if d < k:
                queue.append(u)
                removed.add(u)
        return queue, removed

    @staticmethod
    def prune_nodes(adj, deg_map, queue, removed, k):
        """
        Main peeling loop: pop a node, decrease neighbors' degree,
        and enqueue those that drop below k.
        """
        while queue:
            u = queue.popleft()
            for w in adj[u]:
                if w not in removed:
                    deg_map[w] -= 1
                    if deg_map[w] < k:
                        removed.add(w)
                        queue.append(w)

    @staticmethod
    def collect_remaining_nodes(adj, removed):
        """Return nodes that were NOT pruned away."""
        remaining = set()
        for u in adj:
            if u not in removed:
                remaining.add(u)
        return remaining

    @staticmethod
    def count_components(adj):
        """
        Count how many connected components exist in this adjacency dict.
        """
        visited = set()
        comp_cnt = 0
        for u in adj:
            if u not in visited:
                comp_cnt += 1
                q = deque([u])
                visited.add(u)
                while q:
                    x = q.popleft()
                    for y in adj[x]:
                        if y not in visited:
                            visited.add(y)
                            q.append(y)
        return comp_cnt




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
