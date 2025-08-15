#!/usr/bin/env python3
# Auto-generated for 5514755

STUDENT_ID = "5514755"
STUDENT_NAME = "Xinze Li"

# ======= 学生代码 =======
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
        n = G.vertex_num
        sd = [0] * n


        adj_container = None
        # 1) list-of-lists
        for name, val in G.__dict__.items():
            if isinstance(val, list) and len(val) == n and all(isinstance(x, list) for x in val):
                adj_container = ('list', name)
                break
        # 2) dict-of-lists
        if adj_container is None:
            for name, val in G.__dict__.items():
                if isinstance(val, dict):
                    keys = set(val.keys())
                    if all(isinstance(k, int) for k in keys) and keys.issubset(set(range(n))) \
                       and all(isinstance(val[k], list) for k in keys):
                        adj_container = ('dict', name)
                        break
        # 3) methods
        if adj_container is None:
            if hasattr(G, 'neighbors') and callable(G.neighbors):
                get_nbrs = G.neighbors
            elif hasattr(G, 'get_neighbors') and callable(G.get_neighbors):
                get_nbrs = G.get_neighbors
            elif hasattr(G, 'adjacent') and callable(G.adjacent):
                get_nbrs = G.adjacent
            else:
                raise AttributeError(
                    "'UndirectedUnweightedGraph' has no recognizable adjacency structure"
                )
        else:
            kind, name = adj_container
            if kind == 'list':
                get_nbrs = lambda u, name=name: getattr(G, name)[u]
            else:  # dict
                get_nbrs = lambda u, name=name: getattr(G, name).get(u, [])


        for v in range(n):
            nbrs = get_nbrs(v)
            if len(nbrs) < k:
                continue

            nbr_set = set(nbrs)
            local_deg = {u: 0 for u in nbrs}

            for u in nbrs:
                for w in get_nbrs(u):
                    if w in nbr_set:
                        local_deg[u] += 1


            q = deque(u for u, d in local_deg.items() if d < k)
            while q:
                u = q.popleft()
                if u not in local_deg:
                    continue
                for w in get_nbrs(u):
                    if w in local_deg:
                        local_deg[w] -= 1
                        if local_deg[w] < k:
                            q.append(w)
                del local_deg[u]

            if not local_deg:
                continue


            visited = set()
            cc = 0
            for u in local_deg:
                if u in visited:
                    continue
                cc += 1
                stack = [u]
                visited.add(u)
                while stack:
                    x = stack.pop()
                    for w in get_nbrs(x):
                        if w in local_deg and w not in visited:
                            visited.add(w)
                            stack.append(w)
            sd[v] = cc

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
