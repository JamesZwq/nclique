#!/usr/bin/env python3
# Auto-generated for 5500941

STUDENT_ID = "5500941"
STUDENT_NAME = "Hanchao Xie"

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
        n = G.vertex_num
        tau = [0] * n

        for v in range(n):
            neighbors = set(G.adj_list[v])

            if len(neighbors) == 0:
                tau[v] = 0
                continue

            # Construct adj list for subgraph
            adjNeighbor = {}
            for ni in neighbors:
                adjNeighbor[ni] = []
                for w in G.adj_list[ni]:
                    if w in neighbors:
                        adjNeighbor[ni].append(w)

            tau[v] = kCoreBaseStructuralDiversity.countCores(adjNeighbor, k)

        return tau

    @staticmethod
    def countCores(adjDict, k):
        if not adjDict:
            return 0

        if k == 0:
            return kCoreBaseStructuralDiversity.getCCNum(adjDict)

        cv = kCoreBaseStructuralDiversity.getCoreVertices(adjDict, k)

        if not cv:
            return 0

        adjCore = {}
        for cvi in cv:
            adjCore[cvi] = []
            for v in adjDict.get(cvi, []):
                if v in cv:
                    adjCore[cvi].append(v)

        return kCoreBaseStructuralDiversity.getCCNum(adjCore)

    @staticmethod
    def getCoreVertices(adjDict, k):
        if not adjDict:
            return 0
        
        degrees = {v: len(neighbors) for v, neighbors in adjDict.items()}
        vertices = set(adjDict.keys())

        # Remove vertices with degree less than k
        queue = deque()
        for v in vertices:
            if degrees[v] < k:
                queue.append(v)

        while queue:
            v = queue.popleft()
            if v not in vertices:continue

            vertices.remove(v)

            # Update degrees of neighbor
            for v in adjDict.get(v, []):
                if v in vertices:
                    degrees[v] -= 1
                    if degrees[v] < k:
                        queue.append(v)

        return vertices

    @staticmethod
    # Get Number of connected components
    def getCCNum(adjDict):
        if not adjDict:
            return 0
        
        visited = set()
        components = 0

        for adj in adjDict:
            if adj not in visited:
                # BFS
                components += 1
                queue = deque([adj])
                visited.add(adj)

                while queue:
                    u = queue.popleft()
                    for v in adjDict.get(u, []):
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)

        return components


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
