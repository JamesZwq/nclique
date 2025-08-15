#!/usr/bin/env python3
# Auto-generated for 5537542

STUDENT_ID = "5537542"
STUDENT_NAME = "Jinghua Li"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque,defaultdict
################################################################################

from collections import deque, defaultdict

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

        def get_neighbour_subgraph(v):#Extracts the neighbor-induced subgraph of vertex v
            neighbors=G.adj_list[v]#First, remove all neighbouring nodes of v.
            sub_adj={u:set() for u in neighbors}
            for u in neighbors:#Then iterate through these neighbours u. For each neighbour w of u,if w∈N(v), add edge (u,w).
                for w in G.adj_list[u]:
                    if w in sub_adj:
                        sub_adj[u].add(w)
                        sub_adj[w].add(u)
            return sub_adj#Output a dictionary sub_adj, which is an induced subgraph in adjacency list form.

        def compute_k_core(graph,k):#Computes the k-core of the subgraph using degree peeling.Extract the k-core from the induced subgraph and remove all nodes with degrees less than k.
            degree={u:len(neigh) for u,neigh in graph.items()}#Initialise the degree of each node;
            queue=deque([u for u in graph if degree[u]<k])#Maintain a queue, initially adding all nodes with degrees less than k.
            while queue:#Each time a node u is popped from the queue, remove it from its neighbours.
                u=queue.popleft()
                for v in graph[u]:
                    if v in graph:
                        graph[v].discard(u)
                        if len(graph[v])==k-1:#If neighbour v's degree becomes less than k, it is also added to the queue.
                            queue.append(v)
                del graph[u]#Delete node u until there are no nodes with degree less than k.
            return graph

        def count_components(graph):#Counts connected components in the given adjacency list.
            visited=set()
            count=0
            def bfs(start):
                q=deque([start])
                visited.add(start)
                while q:
                    node=q.popleft()
                    for neighbor in graph[node]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            q.append(neighbor)

            for u in graph:
                if u not in visited:#Iterate through each unvisited node.
                    bfs(u)#For each newly discovered node, perform a BFS once.
                    count+=1#Each time BFS scans an entire connected component, the count is incremented by one.
            return count

        n=G.vertex_num
        τ=[0]*n
        for v in range(n):#main function
            nbr_sub=get_neighbour_subgraph(v)
            k_core_sub=compute_k_core(nbr_sub,k)
            τ[v]=count_components(k_core_sub)
        return τ


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
