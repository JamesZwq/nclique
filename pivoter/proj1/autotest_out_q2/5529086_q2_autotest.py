#!/usr/bin/env python3
# Auto-generated for 5529086

STUDENT_ID = "5529086"
STUDENT_NAME = "Zhigeng Wei"

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
        n = G.vertex_num
        sd = [0] * n
        adj = G.adj_list

        # Get neighbor-induced subgraph's components number
        def count_components(nodes, sub_adj):
            visited = set()
            component_count = 0

            # travel each node
            for node in nodes:
                if node in visited:
                    continue
                # get new component
                component_count += 1
                # BFS by queue
                queue = [node]
                visited.add(node)

                while queue:
                    current = queue.pop(0)
                    # travel each nbrs for this node
                    for neighbor in sub_adj[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)

            return component_count

        # Traverse all v, get all τ_k(v)
        for v in range(n):
            nbrs = adj[v]
            nbrs_num = len(nbrs)
            if nbrs_num < k:
                sd[v] = 0
                continue

            # Step1: Construct the neighbor-induced subgraph: for adj & degree
            nbr_set = set(nbrs)
            # construct sub-graph adj list
            sub_adj = {}
            for u in nbrs:
                sub_adj[u] = []
            # construct sub-graph degree dictionary
            sub_deg = {}
            for u in nbrs:
                deg = 0
                # Make sure that the neighboring points do not include v itself
                commonset = set(adj[u]) & nbr_set
                for w in commonset:
                    sub_adj[u].append(w)
                    deg += 1
                sub_deg[u] = deg
                # print(nbrs)
                # print(sub_adj)
                # print(sub_deg)
                # print()

            # Step2: Handling situations with different k values
            # Situation 1: k == 0
            if k == 0:
                sd[v] = count_components(nbrs, sub_adj)
                continue

            # Situation 2: k > 0
            queue = deque()
            for node in nbrs:
                if sub_deg[node] < k:
                    queue.append(node)
            # Record the removed nodes
            removed_nodes = set()

            # Start the BFS reduce "nodes" process
            while queue:
                current = queue.popleft()
                if current in removed_nodes:
                    continue
                # mark "removed"
                removed_nodes.add(current)

                # Update the degrees of all nbrs
                for nbr in sub_adj[current]:
                    if nbr not in removed_nodes:
                        sub_deg[nbr] = sub_deg[nbr] - 1
                        # Further check if it satisfies "k-core"
                        if sub_deg[nbr] == k - 1:
                            queue.append(nbr)
            # print(queue)
            # print()

            # Step3: Get remain nbrs & remain adj list
            remain_nbrs = []
            for node in nbrs:
                if node not in removed_nodes:
                    remain_nbrs.append(node)
            remain_set = set(remain_nbrs)

            remain_adj = {}
            for node in remain_nbrs:
                remainbrs = []
                for nbr in sub_adj[node]:
                    if nbr in remain_set:
                        remainbrs.append(nbr)
                remain_adj[node] = remainbrs

            # Step4: Caculate sd[v] by components method
            if not remain_nbrs:
                sd[v] = 0
            else:
                sd[v] = count_components(remain_nbrs, remain_adj)

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
