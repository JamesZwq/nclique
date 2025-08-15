#!/usr/bin/env python3
# Auto-generated for 5509423

STUDENT_ID = "5509423"
STUDENT_NAME = "Chi Zhang"

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
        # TODO
        if k == None:
            k = 0
        n = G.vertex_num
        T = [0] * n
        if k < 0 or n == 0:
            return T
        graph = {u: set(G.adj_list[u]) for u in range(n)}
        core = kCoreBaseStructuralDiversity.computeCore(graph)

        availNode = {u for u, c in core.items() if c >= k}
        if not availNode:
            return T
        components = kCoreBaseStructuralDiversity.getComponent(graph, availNode)

        # process the components and its neighbors O（n）
        for component in components:
            related = kCoreBaseStructuralDiversity.findRelated(component, graph)


            kCoreBaseStructuralDiversity.processComponent(
                graph, related, core, k, T
            )

        return T

    @staticmethod
    def findRelated(component, graph):
     # get component itself and neighbors，used for later processing
     # time complexity：O（n)
        if component == None:
            return set()

        related = set(component)
        for u in component:
            related |= graph[u]
        return related





    @staticmethod
    def computeCore(graph):
    # use coreNumber to get the k core
    #    time complexity：O（n+m)
        n = len(graph)
        degree = {u: len(graph[u]) for u in graph}
        if not degree:
            return {}

        maxDegree = max(degree.values())
        num = [0] * (maxDegree + 1)
        for d in degree.values():
            num[d] += 1
        start = [0] * (maxDegree + 1)
        for d in range(1, maxDegree + 1):
            start[d] = start[d - 1] + num[d - 1]

        order = [None] * n
        pos = {}
        for u, d in degree.items():
            idx = start[d]
            order[idx] = u
            pos[u] = idx
            start[d] += 1
        for d in range(maxDegree, 0, -1):
            start[d] = start[d - 1]
        start[0] = 0

        core = kCoreBaseStructuralDiversity.CoreNumber(graph)
        return core
    @staticmethod
    def CoreNumber(graph):
        # count the core layer at each node:
        # count the degrees of each node, put the node with degree
        # refresh the degree of neighbor
        # time complexity：O（n+m)

        degree = {u: len(neighborS) for u, neighborS in graph.items()}
        if not degree:
            return {}
        # find the max degree
        maxDegree = max(degree.values())
        # a list use to save all degree which equal to d
        buckets = [[] for _ in range(maxDegree + 1)]
        for u, d in degree.items():
            if d<0:
                continue
            buckets[d].append(u)

        # from 0 to max degree
        core = {}
        processed = set()
        for tempDegree in range(maxDegree + 1):
            for u in buckets[tempDegree]:
                if u in processed:
                    continue
                processed.add(u)
                core[u] = tempDegree
                for v in graph[u]:
                    if v not in processed and degree[v] > tempDegree:
                        Pre = degree[v]
                        degree[v] -= 1
                        buckets[Pre - 1].append(v)

                # minus the degree if its not visited and degree >tempDegree
                for v in graph[u]:
                    if v not in processed and degree[v] > tempDegree:
                        Pre = degree[v]
                        degree[v] -= 1
                        buckets[Pre - 1].append(v)

        return core

    @staticmethod
    def getComponent(graph, nodes):
        # find component of nodes
        # time complexity：O(n+m)
        visited = set()
        clist = []
        for u in nodes:
            if u in visited:
                continue
            component = {u}
            visited.add(u)
            queue = [u]
            for x in queue:
                for y in graph[x]:
                    # focus on node and never been visited
                    if y in nodes and y not in visited:
                        visited.add(y)
                        component.add(y)
                        queue.append(y)
            clist.append(component)
        return clist


    @staticmethod
    def processComponent(graph, related, core, k, T):
        # find neighbor whose core >=k
        # surbgraph < k , count the component
        # time complexity：O(n)
        availNeighber = {
            v: {u for u in graph[v] if core.get(u, 0) >= k}
            for v in related
        }

        for v in related:
            avaNeighbors = availNeighber[v]
            if len(avaNeighbors) < k:
                continue
            subGragh = {
                u: avaNeighbors & availNeighber[u]
                for u in avaNeighbors
            }

            i = kCoreBaseStructuralDiversity.Kcore(subGragh, k)
            T[v] = max(T[v], i)


    @staticmethod
    def Kcore(subGragh, k):
        # remove all degree < k nodes and return the number of connected remain nodes
       # time complexity：O(n)
        degree = {u: len(subGragh[u]) for u in subGragh}
        queue = [u for u, d in degree.items() if d < k]
        removed = set(queue)

        # use bfs to removed
        for u in queue:
            for w in subGragh[u]:
                if w not in removed:
                    degree[w] -= 1
                    if degree[w] < k:
                        removed.add(w)
                        queue.append(w)

        remianNodes = set(subGragh) - removed
        if not remianNodes:
            return 0
        visited = set()
        componentNum = 0
        for start in remianNodes:
            if start in visited:
                continue
            componentNum += 1
            queue = [start]
            visited.add(start)
            for u in queue:
                for w in subGragh[u]:
                    if w in remianNodes and w not in visited:
                        visited.add(w)
                        queue.append(w)
        return componentNum
# use bucket algorithm to compute each node's core number wit can remove the useless nodes.
# keep only those nodes whose core number >=k, in this surbgraph find the component
# for surbgraph process their neighbors
# finally get the available neighbor and removed the degree < k nodes. and count the remain components.
# the best time complexity: O(n+m)
# the worst time complexity: O(n+m)
# the average time complexity: O(n+m)



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
