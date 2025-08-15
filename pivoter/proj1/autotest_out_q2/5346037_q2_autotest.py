#!/usr/bin/env python3
# Auto-generated for 5346037

STUDENT_ID = "5346037"
STUDENT_NAME = "Elliot Berghofer"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
import copy
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass
    
    
    @staticmethod
    def _neighbour_subgraph(G, curr_v, bitmap): # Note we don't include v in this.
        n = G.vertex_num
        # Shallow copy.
        NSG = copy.copy(G) 
        NSG.adj_list = [[] for _ in range(n)]
        
        # Mark neighbours of v.
        for u in G.adj_list[curr_v]:                                                            
            bitmap[u] = True
        
        # Add all the neighbours for each u where nbr_u \in nbr_v.
        for u in G.adj_list[curr_v]:
            NSG.adj_list[u] = [nbr_u for nbr_u in G.adj_list[u] if bitmap[nbr_u]]
            
        # Unmark neighbours of v.
        for u in G.adj_list[curr_v]:
            bitmap[u] = False
        
        return NSG

    @staticmethod
    def _trackConnectComponents(G, nbr_v, nbr_v_removed):
            visited = nbr_v_removed
            cc = []
            for u in nbr_v:
                if not visited[u]:
                    # DFS (iterative style)
                    stack = deque([u])
                    comp = []
                    while stack:
                        v = stack.pop()
                        if visited[v]:
                            continue
                        visited[v] = True
                        comp.append(v)
                        for x in G.adj_list[v]:
                            stack.append(x)
                    cc.append(comp) 
            return cc

    @staticmethod
    def _kcore(G, k):
        n = G.vertex_num
        
        # 1. Initialise queue with all v with deg(v) < k.
        deg = [len(G.adj_list[v]) for v in range(n)]
        q = deque(v for v, deg_v in enumerate(deg) if deg_v < k)
        
        # 2. Remove all v with deg(v) < k.
        removed = [False] * n
        while q:
            # Pop a vertex from the queue
            v = q.popleft()
            if removed[v]: # v has already been removed.
                continue
            # Remove v from G
            removed[v] = True
            
            # Decrement the degree of each neighbour
            for u in G.adj_list[v]:
                deg[u] -= 1
                # If deg(u) < k push onto the queue.
                if deg[u] < k:
                    q.append(u)
        
        return removed

    @staticmethod
    def _kcore_components(G, k, nbr_v):
        # 1. Run k-core.
        removed = kCoreBaseStructuralDiversity._kcore(G, k)
        # 2. Return the Connected Components.
        kcores = kCoreBaseStructuralDiversity._trackConnectComponents(G, nbr_v, removed)
        return kcores
    
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
        
        # 1. Global k-core, remove all nodes not in the k-core O(n + m).
        removed = kCoreBaseStructuralDiversity._kcore(G, k)
        #    Copy the smaller graph. O(n + m)
        newG = copy.copy(G)
        newG.adj_list = [[] for _ in range(n)]
        for v in range(n):
            if not removed[v]:
                for u in G.adj_list[v]:
                    if not removed[u]:
                        newG.adj_list[v].append(u)
        G = newG
        #    Only nodes that are in k-cores exist now.
        
        bitmap = [False] * n
        for v in range(n):
            # For n nodes max degree is (n-1) so k-core doesn't exist if (n-1) < k.
            if len(G.adj_list[v])-1 >= k:
                NSG   = kCoreBaseStructuralDiversity._neighbour_subgraph(G, v, bitmap)             # Find Neighbour-induced Subgraph.
                sd[v] = len(kCoreBaseStructuralDiversity._kcore_components(NSG, k, G.adj_list[v])) # Compute k-core-based structural diversity for node v.

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
