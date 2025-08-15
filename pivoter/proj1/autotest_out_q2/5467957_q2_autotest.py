#!/usr/bin/env python3
# Auto-generated for 5467957

STUDENT_ID = "5467957"
STUDENT_NAME = "Yaoyu Xie"

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
        result = []
        for v in range(G.vertex_num):
            neighbors = set(G.adj_list[v])
            if not neighbors:
                result.append(0)
                continue
            induced_adj = kCoreBaseStructuralDiversity._induce_adjacent_map(G, neighbors)
            diversity = kCoreBaseStructuralDiversity._find_kcore_components(induced_adj, k)
            result.append(diversity)
        return result

    @staticmethod
    def _induce_adjacent_map(G, node_set):
        """从 G 中抽取节点集合 node_set 的邻接图"""
        return {
            u: [v for v in G.adj_list[u] if v in node_set]
            for u in node_set
        }

    @staticmethod
    def _find_kcore_components(adj, k):
        if not adj:
            return 0

        # 初始化每个节点的度
        deg_map = {u: len(v_list) for u, v_list in adj.items()}
        to_remove = set(u for u, deg in deg_map.items() if deg < k)

        # 剥除所有不满足 k-core 条件的节点
        queue = deque(to_remove)
        while queue:
            curr = queue.popleft()
            for nbr in adj[curr]:
                if nbr not in to_remove:
                    deg_map[nbr] -= 1
                    if deg_map[nbr] < k:
                        to_remove.add(nbr)
                        queue.append(nbr)

        # 统计剩余图中的连通块数量
        valid_nodes = set(adj) - to_remove
        if not valid_nodes:
            return 0

        visited = set()
        kcore_count = 0

        while valid_nodes:
            node = valid_nodes.pop()
            kcore_count += 1
            stack = [node]
            visited.add(node)
            while stack:
                curr = stack.pop()
                for neighbor in adj[curr]:
                    if neighbor not in visited and neighbor in valid_nodes:
                        visited.add(neighbor)
                        stack.append(neighbor)
                        valid_nodes.discard(neighbor)

        return kcore_count

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
