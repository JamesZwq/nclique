#!/usr/bin/env python3
# Auto-generated for 5473354

STUDENT_ID = "5473354"
STUDENT_NAME = "Yuzhe Xia"

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
        
        # Build a subgraph to calculate k-cores
        def subgraph(G, node, neighbor_nodes):
            graph = {}
            for neighbor_node in neighbor_nodes:
                graph[neighbor_node] = [node for node in G.adj_list[neighbor_node] if node in neighbor_nodes]
            return graph
        
        def calc_kcores(G, k, neighbor_nodes):
            degree = {node: len(G[node]) for node in G}
            remove_nodes = [node for node in neighbor_nodes if degree[node] < k]
            # Remove nodes which degree < k
            while remove_nodes:
                remove_node = remove_nodes.pop(0)
                degree[remove_node] = 0
                for neighbor in G[remove_node]:
                    if degree[neighbor] > 0:
                        degree[neighbor] -= 1
                        if degree[neighbor] == k-1:
                            degree[neighbor] = 0
                            remove_nodes.append(neighbor)
            # Get the rest nodes and find the connected components
            node_rest = [node for node in neighbor_nodes if degree[node] >=k]
            visited = set()
            core_count = 0
            for node in node_rest:
                if node in visited:
                    continue
                q = deque()
                q.append(node)
                # New connected component
                core_count += 1
                while q:
                    now_node = q.popleft()
                    if now_node in visited:
                        continue
                    visited.add(now_node)
                    for neighbor in G[now_node]:
                        if degree[neighbor] == 0 or neighbor in visited:
                            continue
                        q.append(neighbor)
            
            return core_count
        
        n = G.vertex_num
        sd = [0] * n
        for i in range(G.vertex_num):
            neighbor_nodes = set(G.adj_list[i])
            if not neighbor_nodes:
                sd[i] = 0
            else:
                G_sub = subgraph(G, i, neighbor_nodes)
                sd[i] = calc_kcores(G_sub, k, neighbor_nodes)
        return sd

    
    # Time complexity: O(|V| * (|V|+|E|))
    # Space complexity: O(|V|+|E|) to store the graph
    
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
