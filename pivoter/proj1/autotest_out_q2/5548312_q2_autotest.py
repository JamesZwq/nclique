#!/usr/bin/env python3
# Auto-generated for 5548312

STUDENT_ID = "5548312"
STUDENT_NAME = "Jinchuan He"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    def k_cores(self, G):
        """
        Use flat array calculate the k-core number of each vertexes in G.
        Returns
        -------
        d[v]: is the k-core number of v.
        """
        n = G.vertex_num
        if n == 0:
            return []
        
        # 1. init degrees
        d = [len(adj) for adj in G.adj_list]
        core = d.copy()
        
        # 2. sort the nodes by degree
        nodes = sorted(range(n), key=lambda x: (d[x], x)) # D
        pos = {u: i for i, u in enumerate(nodes)}  # p
        
        # 3. init bin_boundaries[d] = degree < d +1
        max_degree = max(d) if d else 0
        bin_boundaries = [0] * (max_degree + 1) # b
        for deg in range(1, max_degree + 1):
            # find the first degree >= deg node index
            bin_boundaries[deg] = bin_boundaries[deg - 1]
            while bin_boundaries[deg] < n and d[nodes[bin_boundaries[deg]]] < deg:
                bin_boundaries[deg] += 1
        
        # 4. calculate cores of k
        for i in range(n):
            v = nodes[i]
            for u in G.adj_list[v]:
                if core[u] > core[v]:
                    du = core[u]
                    pos_u = pos[u]
                    if pos_u >= bin_boundaries[du]:
                        w = nodes[bin_boundaries[du]]
                        if u != w:
                            pos[u], pos[w] = bin_boundaries[du], pos_u
                            nodes[pos_u], nodes[bin_boundaries[du]] = nodes[bin_boundaries[du]], nodes[pos_u]
                        bin_boundaries[du] += 1
                        core[u] -= 1
        return core

    def subgraph(self, G, nodes: list) -> list:
        """
        generate subgraph of G, contains only nodes
        """

        limits = set(nodes)
        nodes.sort()

        nodes_map = {}
        for i in range(len(nodes)):
            v = nodes[i]
            nodes_map[v] = i
        # print("nodes map:", nodes_map)

        edges = []
        for i in range(len(nodes)):
            v = nodes[i]
            for u in G.adj_list[v]:
                if u not in limits:
                    continue

                # rename the u to j
                j = nodes_map[u]

                if i == j:
                    continue

                if i < j:
                    edges.append((i, j))
                else:
                    edges.append((j, i))
        edges = list(set(edges))
        edges.sort()
        
        # print("sub graphs", edges)

        return [(len(nodes), len(edges))] + edges


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

        # init parameters
        n = G.vertex_num
        kcsd = kCoreBaseStructuralDiversity()

        # 1. calculate all nbr for G
        nbrs = [[] for _ in range(n)]
        for v in range(n):
            # read all neighbors of v, 
            # and select the edges that contain these points
            nvs = list(set(G.adj_list[v]))
            if len(nvs) > 0:
                nbrs[v] = kcsd.subgraph(G, nvs)

        # 2. calculate K(nbr,k=k), and count the tau(k=k), save to sd
        sd = [0] * n
        for v in range(n):

            if len(nbrs[v]) == 0:
                continue

            G1 = UndirectedUnweightedGraph(nbrs[v])
            k_core = kcsd.k_cores(G1)

            # filter the degree >= k vertexes
            uList = []
            for u, degree in enumerate(k_core):
                if degree >= k:
                    uList.append(u)

            # use BFS get all vertexes in current group
            cores = []
            if len(uList) > 0:
                stack = uList.copy()
                while len(stack) > 0:
                    u = stack.pop()
                    group = []
                    queue = []
                    queue.append(u)
                    while len(queue) > 0:
                        curU = queue.pop()
                        group.append(curU)
                        neighbors = G1.adj_list[curU]
                        neighbors = set(neighbors).intersection(set(uList)) # filter deleted
                        neighbors = list(filter(lambda x: x not in group and x not in queue, neighbors))
                        queue = queue + neighbors

                    cores.append(list(set(group)))
                    stack = list(set(stack).difference(set(group)))

            sd[v] = len(cores)

        return sd


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
