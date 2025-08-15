#!/usr/bin/env python3
# Auto-generated for 5597260

STUDENT_ID = "5597260"
STUDENT_NAME = "Jiangkuan Yu"

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
        

        fulldegrees = [0] * n
        degrees = [len(G.adj_list[u]) for u in range(n)]
        
        queue = deque()
        removed = [False] * n
        
        for u in range(n):
            if degrees[u] == 0:
                queue.append(u)
                removed[u] = True
        
        current_core = 0
        while current_core < n:
            while queue:
                u = queue.popleft()
                fulldegrees[u] = current_core
                for neighbor in G.adj_list[u]:
                    if not removed[neighbor]:
                        degrees[neighbor] -= 1
                        if degrees[neighbor] <= current_core:
                            queue.append(neighbor)
                            removed[neighbor] = True
            
            current_core += 1
            for u in range(n):
                if not removed[u] and degrees[u] <= current_core:
                    queue.append(u)
                    removed[u] = True
            
            if not queue:
                for u in range(n):
                    if not removed[u]:
                        fulldegrees[u] = degrees[u]
                break
        
        def find_k_nbr(neighbors):
            if len(neighbors) < k:
                return 0
            
            valid_neighbors = set(u for u in neighbors if fulldegrees[u] >= k)
            if len(valid_neighbors) < k:
                return 0
            
            degrees = {}
            for u in valid_neighbors:
                degrees[u] = 0
                for w in G.adj_list[u]:
                    if w in valid_neighbors:
                        degrees[u] += 1
            
            queue = deque()
            removed = set()
            
            for u in valid_neighbors:
                if degrees[u] < k:
                    queue.append(u)
                    removed.add(u)
            
            while queue:
                u = queue.popleft()
                for w in G.adj_list[u]:
                    if w in valid_neighbors and w not in removed:
                        degrees[w] -= 1
                        if degrees[w] < k:
                            queue.append(w)
                            removed.add(w)
            
            k_core_nodes = valid_neighbors - removed
            
            if not k_core_nodes:
                return 0
            
            visited = set()
            components = 0
            k_core_set = set(k_core_nodes)
            
            def bfs(node):
                visited.add(node)
                queue = deque([node])
                while queue:
                    node = queue.popleft()
                    for neighbor in G.adj_list[node]:
                        if neighbor in k_core_set and neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
            
            for node in k_core_nodes:
                if node not in visited:
                    bfs(node)
                    components += 1
            
            return components
        
        for v in range(n):
            neighbors = set(G.adj_list[v])
            sd[v] = find_k_nbr(neighbors)
        
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
