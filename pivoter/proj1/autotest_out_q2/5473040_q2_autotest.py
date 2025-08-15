#!/usr/bin/env python3
# Auto-generated for 5473040

STUDENT_ID = "5473040"
STUDENT_NAME = "Jordan Beard"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    def componentDiscovery(self, G, k, v):
        """
        Parameters
        ----------
        G : UndirectedUnweightedGraph
        k : int
        v : int
        Returns
        -------
        int  # τ_k(v) for given v
        """

        #if no neighbours then count is 0, O(n)
        if G.adj_list[v] == [] or k > G.vertex_num or len(G.adj_list[v]) < k:
            return 0

        #O(1) so for overall loop O(n)
        n = G.vertex_num

        #get neighbours of vertex being checked, O(n) for a connected graph
        #set to check membership in O(1) so
        neighbourhood = {u for u in G.adj_list[v]}

        #O(n) in worst case so for all n O(n^2)
        if len(neighbourhood) == 0:
            return 0

        #hash map to access and count vertex count in O(1)
        vertexDeg = defaultdict(int)

        #for all vertices connected to vertex being checked, v, count degree when considering only neighbours of v
        #neighbourhood size n in fully connected graph so O(n)*O(n) = O(n^2) over all v ~O(n^3)
        for vertex in neighbourhood:
            vertexDeg[vertex] = len([u for u in G.adj_list[vertex] if u in neighbourhood])

        #create q of vertices with degree less than k, O(n) in the worst case
        q = deque([vertex, degree] for vertex, degree in vertexDeg.items() if degree < k)

        #while loop O(n) if all vertices are less than k
        while(q):

            #pop and double check eligible neighbour
            u = q.popleft()
            if u[0] not in neighbourhood:
                continue

            #remove vertex from valid graph
            neighbourhood.remove(u[0])

            #decrement neighbour count for neighbours, if any drop below k add to q
            #if fully connected O(n) -> O(n^2) in while loop
            for w in G.adj_list[u[0]]:

                if w in vertexDeg.keys() and w in neighbourhood:
                    vertexDeg[w] -= 1

                    if vertexDeg[w] < k:
                        q.append((w, vertexDeg[w]))

        #if neighbours of vertex are gone then 0
        if len(neighbourhood) == 0:
            return 0

        #Make sure vertex being checked not in neighbourhood
        #O(1)
        if v in neighbourhood:
            neighbourhood.remove(v)

        #DFS to count components
        #algorithm taken from GeeksForGeeks
        visited = [False] * n
        count = 0

        #component count BFS O(n + m)
        for vertex in neighbourhood:
            if not visited[vertex]:
                count += 1
                self.DFS(vertex, visited, neighbourhood, G)

        return count

    def DFS(self, vertex, visited, neighbourhood, G):
        visited[vertex] = True
        for w in G.adj_list[vertex]:
            if not visited[w] and w in neighbourhood:
                self.DFS(w, visited, neighbourhood, G)


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

        SD = kCoreBaseStructuralDiversity()

        #n vertices, componentDiscovery O(n^2) so total O(n^3)
        for v in range(n):
            sd[v] = SD.componentDiscovery(G, k, v)

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
