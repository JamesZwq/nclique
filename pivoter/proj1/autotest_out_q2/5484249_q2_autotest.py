#!/usr/bin/env python3
# Auto-generated for 5484249

STUDENT_ID = "5484249"
STUDENT_NAME = "Ruiman Xu"

# ======= 学生代码 =======
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
        List[int]  # τ_k(v) for all vertices v
        """
        def build_induced_subgraph(neighbors):
            # Build the neighbor-induced subgraph from N(v)
            subgraph = {}
            for u in neighbors:
                subgraph[u] = set()
                for v in G.adj_list[u]:
                    if v in neighbors:
                        subgraph[u].add(v)
            return subgraph

        def compute_k_core(subgraph, k):
            # Iteratively remove nodes with degree < k
            subgraph = {u: set(adj) for u, adj in subgraph.items()}  # deep copy
            degrees = {u: len(adj) for u, adj in subgraph.items()}
            changed = True
            while changed:
                changed = False
                remove = [u for u in subgraph if degrees[u] < k]
                if remove:
                    changed = True
                    for u in remove:
                        for v in subgraph[u]:
                            subgraph[v].discard(u)
                            degrees[v] -= 1
                        del subgraph[u]
                        del degrees[u]
            return subgraph

        def count_connected_components(subgraph):
            # Count connected components using DFS
            visited = set()
            count = 0

            def dfs(u):
                stack = [u]
                while stack:
                    node = stack.pop()
                    if node not in visited:
                        visited.add(node)
                        stack.extend(subgraph[node] - visited)

            for u in subgraph:
                if u not in visited:
                    dfs(u)
                    count += 1
            return count

        # Main logic: compute τ_k(v) for each vertex v
        τ = [0] * G.vertex_num
        for v in range(G.vertex_num):
            neighbors = set(G.adj_list[v])
            if not neighbors:
                τ[v] = 0
                continue

            subgraph = build_induced_subgraph(neighbors)
            k_core_subgraph = compute_k_core(subgraph, k)
            τ[v] = count_connected_components(k_core_subgraph)

        return τ

    @staticmethod
    def explain_k_core_diversity(G, v, max_k=4):
        """
        Print detailed explanation of τ_k(v) for k = 0 to max_k
        """
        def build_induced_subgraph(neighbors):
            # Build the neighbor-induced subgraph
            subgraph = {}
            for u in neighbors:
                subgraph[u] = set()
                for x in G.adj_list[u]:
                    if x in neighbors:
                        subgraph[u].add(x)
            return subgraph

        def compute_k_core(subgraph, k):
            # Prune subgraph to obtain k-core
            subgraph = {u: set(adj) for u, adj in subgraph.items()}
            degrees = {u: len(adj) for u, adj in subgraph.items()}
            changed = True
            while changed:
                changed = False
                remove = [u for u in subgraph if degrees[u] < k]
                if remove:
                    changed = True
                    for u in remove:
                        for v in subgraph[u]:
                            subgraph[v].discard(u)
                            degrees[v] -= 1
                        del subgraph[u]
                        del degrees[u]
            return subgraph

        def get_connected_components(subgraph):
            # Return connected components in the subgraph
            visited = set()
            components = []

            def dfs(u, comp):
                stack = [u]
                while stack:
                    node = stack.pop()
                    if node not in visited:
                        visited.add(node)
                        comp.append(node)
                        stack.extend(subgraph[node] - visited)

            for u in subgraph:
                if u not in visited:
                    comp = []
                    dfs(u, comp)
                    components.append(sorted(comp))
            return components

        neighbors = set(G.adj_list[v])
        print(f"N({v}) = {sorted(neighbors)}\n")

        for k in range(0, max_k + 1):
            subgraph = build_induced_subgraph(neighbors)
            k_core_sub = compute_k_core(subgraph, k)
            components = get_connected_components(k_core_sub)

            print(f"k = {k}")
            print(f"τ_{k}({v}) = {len(components)}")
            if components:
                print(f"𝒦_{k}(nbr_{v}) = {components}\n")
            else:
                print(f"𝒦_{k}(nbr_{v}) = ∅ (no {k}-core)\n")

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
