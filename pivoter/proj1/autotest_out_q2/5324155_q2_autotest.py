#!/usr/bin/env python3
# Auto-generated for 5324155

STUDENT_ID = "5324155"
STUDENT_NAME = "Chenkai Shen"

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
        sd = [0] * n
        
        for v in range(n):
            neighbors = G.adj_list[v]
            
            if len(neighbors) == 0:
                sd[v] = 0
                continue
            
            if len(neighbors) < k + 1:
                sd[v] = 0
                continue
            
            # construct the neighbor subgraph
            node_to_index = {node: i for i, node in enumerate(neighbors)}
            neighbors_degrees = [0] * len(neighbors)
            neighbors_set = set(neighbors)

            adj_list = [[] for _ in range(len(neighbors))]

            for u in neighbors:
                for w in G.adj_list[u]:
                    if w in neighbors_set and u < w:
                        u_idx = node_to_index[u]
                        w_idx = node_to_index[w]
                        adj_list[u_idx].append(w_idx)
                        adj_list[w_idx].append(u_idx)
                        neighbors_degrees[u_idx] += 1
                        neighbors_degrees[w_idx] += 1
            
            # find the k-core of the neighbor subgraph
            queue = deque()
            is_removed = [False] * len(neighbors)
            # initialize the queue
            for i in range(len(neighbors)):
                if neighbors_degrees[i] < k:
                    queue.append(i)
                    is_removed[i] = True
            
            # pelling the graph
            while len(queue) > 0:
                u_idx = queue.popleft()
                for v_idx in adj_list[u_idx]:
                    if not is_removed[v_idx]:
                        neighbors_degrees[v_idx] -= 1
                        if neighbors_degrees[v_idx] < k:
                            queue.append(v_idx)
                            is_removed[v_idx] = True
            
            # count the number of connected components in the graph
            num_connected_components = 0
            visited = [False] * len(neighbors)

            for i in range(len(neighbors)):
                if not visited[i] and not is_removed[i]:
                    num_connected_components += 1

                    temp_queue = deque()
                    temp_queue.append(i)
                    visited[i] = True

                    while len(temp_queue) > 0:
                        curr_idx = temp_queue.popleft()
                        for v_idx in adj_list[curr_idx]:
                            if not visited[v_idx] and not is_removed[v_idx]:
                                temp_queue.append(v_idx)
                                visited[v_idx] = True
            sd[v] = num_connected_components
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
