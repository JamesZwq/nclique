#!/usr/bin/env python3
# Auto-generated for 5527937

STUDENT_ID = "5527937"
STUDENT_NAME = "Guo Li"

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
            输入的无向无权图，含属性
              - G.vertex_num: 顶点总数（0..vertex_num-1）
              - G.adj_list: 长度为 vertex_num 的邻接表列表
        k : int
            指定的 k 值

        Returns
        -------
        List[int]
            长度为 G.vertex_num 的列表，t[v] = τ_k(v)
        """
        n = G.vertex_num
        adj = G.adj_list
        t = [0] * n

        for v in range(n):
            neigh = adj[v]
            if not neigh:
                t[v] = 0
                continue

            # 1) 构造邻居诱导子图的邻接表
            N_set = set(neigh)
            H_adj = {u: [] for u in neigh}
            for u in neigh:
                # 只保留指向同一邻居集内的边
                for w in adj[u]:
                    if w in N_set:
                        H_adj[u].append(w)

            # 2) k-核剥离
            deg = {u: len(H_adj[u]) for u in neigh}
            queue = deque(u for u, d in deg.items() if d < k)
            while queue:
                u = queue.popleft()
                if deg[u] < 0:
                    # 已剥离过
                    continue
                # 从每个剩余邻居中移除 u
                for w in H_adj[u]:
                    if deg[w] >= k:
                        deg[w] -= 1
                        if deg[w] == k - 1:
                            queue.append(w)
                # 标记 u 已移除
                deg[u] = -1

            # 3) 收集剩余节点
            survivors = [u for u, d in deg.items() if d >= k]
            if not survivors:
                t[v] = 0
                continue

            # 4) 统计剩余子图的连通分量数
            seen = set()
            comps = 0
            for u in survivors:
                if u in seen:
                    continue
                comps += 1
                # BFS 标记此分量所有顶点
                dq = deque([u])
                seen.add(u)
                while dq:
                    x = dq.popleft()
                    for w in H_adj[x]:
                        # 仅在 survivors 且未访问时扩展
                        if w not in seen and deg.get(w, -1) >= k:
                            seen.add(w)
                            dq.append(w)
            t[v] = comps

        return t

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
