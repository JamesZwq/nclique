#!/usr/bin/env python3
# Auto-generated for 5513691

STUDENT_ID = "5513691"
STUDENT_NAME = "Lee Zhao"

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

        # The number of nodes in the graph
        n = G.vertex_num

        # Create the set representation of the neighbours of each node
        neighbourSets = [set(adj) for adj in G.adj_list]

        # A list to store the structural diversity value of each node
        diversities = [0] *n

        for node in range(n):

            # Extract the neighbour set of this node, which is the set of all nodes in its ego net
            neighbourSet = neighbourSets[node]

            # If it is an isolated node, the corresponding structural diversity value is 0
            if len(neighbourSet) == 0:
                diversities[node] = 0
                continue

            # Calculate node degrees in ego net
            egoDegrees = {}
            for neighbour in neighbourSet:
                egoDegrees[neighbour] = len(neighbourSets[neighbour] & neighbourSet)

            # K-Core decomposition
            deletedSet = set()
            queue = deque()

            # Collect nodes with ego degress less than k
            for neighbour in neighbourSet:
                if egoDegrees[neighbour] < k:
                    deletedSet.add(neighbour)
                    queue.append(neighbour)

            # Iteratively delete nodes with ego degrees less than k
            while queue:
                currentNode = queue.popleft()

                for neighbour in neighbourSets[currentNode]:

                    # The ego degree of ego neighbours of current node is reduced by 1
                    # Only process this neighbour if it is in the ego net and will not be deleted
                    if neighbour in neighbourSet and neighbour not in deletedSet:
                        egoDegrees[neighbour] -= 1

                        # Check the new ego degree of this neighbour is less than k or not
                        if egoDegrees[neighbour] < k:
                            # If so, we need to delete this node as well
                            deletedSet.add(neighbour)
                            queue.append(neighbour)

            # After deleting nodes according to the k-ego-degree threshold, the remaining nodes are candidates of k-cores
            remainingSet = neighbourSet - deletedSet

            # Calculate the number of k-cores, which is the number of connected componenets in the net composed by the remaining nodes
            if not remainingSet:
                diversities[node] = 0
            else:
                visitedSet = set()
                componentNum = 0

                # Conduct BFS for connected component detection
                for remainingNode in remainingSet:
                    if remainingNode not in visitedSet:
                        componentNum += 1
                        queue = deque([remainingNode])
                        visitedSet.add(remainingNode)

                        while queue:
                            currentNode = queue.popleft()
                            for neighbour in neighbourSets[currentNode]:
                                if neighbour in remainingSet and neighbour not in visitedSet:
                                    visitedSet.add(neighbour)
                                    queue.append(neighbour)

                diversities[node] = componentNum

        return diversities


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
