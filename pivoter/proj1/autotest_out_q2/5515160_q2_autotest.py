#!/usr/bin/env python3
# Auto-generated for 5515160

STUDENT_ID = "5515160"
STUDENT_NAME = "Tianshi Song"

# ======= 学生代码 =======
# ################################################################################
# # You can import any Python Standard Library modules~
# from collections import deque
# ################################################################################

# class kCoreBaseStructuralDiversity(object):
#     def __init__(self):
#         pass

#     staticmethod
#     def process(G, k):
#         """
#         Parameters
#         ----------
#         G : UndirectedUnweightedGraph
#         k : int
#         Returns
#         -------
#         List[int]  # τ_k(v) for all v
#         """
#         # TODO
#         n = G.vertex_num
#         sd = [0] * n
#         return sd


#     ################################################################################
#     # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
#     ################################################################################

from collections import deque

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
            sd[v] = kCoreBaseStructuralDiversity._get_kcore_count(G, v, k)

        return sd

    @staticmethod
    def _get_kcore_count(G, v, k):
        neighbours = G.adj_list[v]
        if not neighbours or len(neighbours) < k:
            return 0

        sub_G = kCoreBaseStructuralDiversity._build_induced_subgraph(G, neighbours)
        return kCoreBaseStructuralDiversity._count_k_core(sub_G, k)

    @staticmethod
    def _build_induced_subgraph(G, neighbour_list):
        node_set = set(neighbour_list)
        sub_G = {u: [] for u in node_set}

        for u in node_set:
            sub_G[u] = [w for w in G.adj_list[u] if w in node_set]

        return sub_G

    @staticmethod
    def _count_k_core(sub_G, k):
        if not sub_G:
            return 0

        degrees = {u: len(neigh) for u, neigh in sub_G.items()}
        deleted = set()
        queue = kCoreBaseStructuralDiversity._initialize_queue(sub_G, degrees, deleted, k)
        kCoreBaseStructuralDiversity._iterative_prune(sub_G, degrees, deleted, queue, k)

        remains = [u for u in sub_G if u not in deleted]
        if not remains:
            return 0

        return kCoreBaseStructuralDiversity._count_components(sub_G, remains, deleted)

    @staticmethod
    def _initialize_queue(sub_G, degrees, deleted, k):
        q = deque([u for u in sub_G if degrees[u] < k])
        deleted.update(q)
        return q

    @staticmethod
    def _iterative_prune(sub_G, degrees, deleted, queue, k):
        while queue:
            u = queue.popleft()
            for v in sub_G[u]:
                if v in deleted:
                    continue
                degrees[v] -= 1
                if degrees[v] < k:
                    deleted.add(v)
                    queue.append(v)

    @staticmethod
    def _count_components(sub_G, remains, deleted):
        visited = set()
        count = 0
        for u in remains:
            if u in visited:
                continue
            count += 1
            kCoreBaseStructuralDiversity._bfs_component(sub_G, u, visited, deleted)
        return count

    @staticmethod
    def _bfs_component(sub_G, start, visited, deleted):
        queue = deque([start])
        visited.add(start)
        while queue:
            curr = queue.popleft()
            for v in sub_G[curr]:
                if v not in deleted and v not in visited:
                    visited.add(v)
                    queue.append(v)

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
