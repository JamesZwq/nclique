#!/usr/bin/env python3
# Auto-generated for 5454598

STUDENT_ID = "5454598"
STUDENT_NAME = "Xiangyi Wang"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
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

        for v in range(n):
            neighbours = G.adj_list[v]
            if not neighbours:
                continue

            neighbours_set = set(neighbours)


            sub_G = {
                u: [w for w in G.adj_list[u] if w in neighbours_set]
                for u in neighbours_set
            }


            sd[v] = kCoreBaseStructuralDiversity._compute_tau_k(sub_G, k)

        return sd

    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def _compute_tau_k(sub_G, k):
        def _prune_k_core(graph, k):
            degrees = {u: len(neigh) for u, neigh in graph.items()}
            removed = set()
            queue = [u for u in graph if degrees[u] < k]
            removed.update(queue)
            front = 0

            while front < len(queue):
                u = queue[front]
                front += 1
                for v in graph[u]:
                    if v not in removed:
                        degrees[v] -= 1
                        if degrees[v] < k:
                            queue.append(v)
                            removed.add(v)


            new_graph = {}
            for u in graph:
                if u not in removed:
                    new_graph[u] = []
                    for v in graph[u]:
                        if v not in removed:
                            new_graph[u].append(v)
            return new_graph

        def _count_components(graph):
            visited = set()
            count = 0
            for node in graph:
                if node not in visited:
                    count += 1
                    queue = [node]
                    front = 0
                    visited.add(node)
                    while front < len(queue):
                        curr = queue[front]
                        front += 1
                        for nei in graph[curr]:
                            if nei not in visited:
                                visited.add(nei)
                                queue.append(nei)
            return count

        if not sub_G:
            return 0

        pruned_graph = _prune_k_core(sub_G, k)
        return _count_components(pruned_graph)

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
