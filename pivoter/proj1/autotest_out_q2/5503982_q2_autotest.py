#!/usr/bin/env python3
# Auto-generated for 5503982

STUDENT_ID = "5503982"
STUDENT_NAME = "Varshith Yadugani"

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
        num_vertices = G.vertex_num
        #print(num_vertices)
        adj = G.adj_list
        #print(adj)
        tau = [0] * num_vertices
        #print(tau)

        # I reuse these arrays for every v to avoid reallocating
        in_set   = [False] * num_vertices
        gone     = [False] * num_vertices
        seen     = [False] * num_vertices
        deg_local= [0] * num_vertices
        touched  = []


        # I count connected components among the still‑active neighbors
        def count_comps(neigh):
            comps = 0
            for start in neigh:
              #print(start)
                if gone[start] or seen[start]:
                    continue
                comps += 1
                #print(comps)
                stack = [start]
                seen[start] = True
                while stack:
                    cur = stack.pop()
                    for nxt in adj[cur]:
                        if not in_set[nxt]:
                            continue
                        if gone[nxt]:
                            continue
                        if seen[nxt]:
                            continue
                        seen[nxt] = True
                        stack.append(nxt)
            return comps

        for v in range(num_vertices):
            neigh = adj[v]
            d = len(neigh)
            #print(neigh)
            #print(d)
            # I skip vertices with no neighbors
            if d == 0:
                continue
            # If k>0 and the neighbor set is too small
            if k > 0 and d <= k:
                continue

            # I mark the neighbor set and reset flags
            for u in neigh:
                in_set[u] = True
                #print(in_set[u])
                gone[u] = False
                seen[u] = False
                touched.append(u)
                #print(touched,end=" ")

            if k == 0:
                # No peeling needed, just count components
                tau[v] = count_comps(neigh)
            else:
                # I compute internal degrees inside N(v)
                for u in neigh:
                  #print("u is:",u)
                    cnt = 0
                    for w in adj[u]:
                      #print("w is:",w)
                        if in_set[w]:
                            cnt += 1
                    deg_local[u] = cnt

                # I peel out vertices whose degree < k
                q = deque()
                for u in neigh:
                  #print(u)
                    if deg_local[u] < k:
                        q.append(u)
                #print(q,end=" ")

                while q:
                    u = q.popleft()
                    #print(u)
                    if gone[u]:
                        continue
                    gone[u] = True
                    for w in adj[u]:
                      #print(w)
                        if in_set[w] and not gone[w]:
                            deg_local[w] -= 1
                            if deg_local[w] == k - 1:
                                q.append(w)
                    #print(q,end=" ")

                # Now I count how many components remain
                tau[v] = count_comps(neigh)

            # I unmark everyone I touched so I can reuse the arrays
            for u in touched:
              #print(u)
              #print("touch is:",touched)
                in_set[u] = False
            touched.clear()

        return tau


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
