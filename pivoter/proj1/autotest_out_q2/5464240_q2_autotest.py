#!/usr/bin/env python3
# Auto-generated for 5464240

STUDENT_ID = "5464240"
STUDENT_NAME = "Taiyu Du"

# ======= 学生代码 =======
class kCoreBaseStructuralDiversity(object):
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

        def find_components(neigh_list):
            """Find all connected groups (components) in the neighbor subgraph"""
            not_visited = set(neigh_list)
            groups = []

            while not_visited:
                start = not_visited.pop()
                stack = [start]
                group = [start]

                while stack:
                    now = stack.pop()
                    for friend in G.adj[now]:
                        if friend in not_visited:
                            stack.append(friend)
                            group.append(friend)
                            not_visited.remove(friend)

                groups.append(group)
            return groups

        def check_k_core(group):
            """Check if all nodes in the group have degree >= k (only inside the group)"""
            deg = {node: 0 for node in group}
            for node in group:
                for friend in G.adj[node]:
                    if friend in group:
                        deg[node] += 1

            print(f"  Node degrees: {deg}")
            return all(deg[node] >= k for node in group)

        result = []
        for v in range(G.vertex_num):
            neigh = G.adj[v]
            neigh_list = list(neigh)
            print(f"\nVertex {v} has neighbors: {neigh_list}")

            groups = find_components(neigh_list)
            print(f"  Found {len(groups)} group(s): {groups}")

            count = 0
            for i, group in enumerate(groups):
                print(f"    Checking group {i}: {group}")
                if check_k_core(group):
                    print(f"  😍 Group {i} is a {k}-core")
                    count += 1
                else:
                    print(f"  🥲 Group {i} is NOT a {k}-core")

            result.append(count)
            print(f"  τ_{k}({v}) = {count}")
        return result

# Run the function and print results
τ = kCoreBaseStructuralDiversity.process(G, 2)
print("Final τ_k values:", τ)

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
