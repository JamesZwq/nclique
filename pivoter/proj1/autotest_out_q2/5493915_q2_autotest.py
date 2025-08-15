#!/usr/bin/env python3
# Auto-generated for 5493915

STUDENT_ID = "5493915"
STUDENT_NAME = "Boqian Hu"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def zero_list(n):
        return [0] * n

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
        tau = [0] * n
        if n == 0:
            return tau

        adj = G.adj_list

        # vertex in current neighborhood
        mark_inN = kCoreBaseStructuralDiversity.zero_list(n)
        # removed by peeling
        removed = kCoreBaseStructuralDiversity.zero_list(n)
        # visited when counting components
        seen = kCoreBaseStructuralDiversity.zero_list(n)
        # degree inside neighborhood
        deg_local = kCoreBaseStructuralDiversity.zero_list(n)
        # increasing timestamp
        cur_tag = 1
        for v in range(n):
            N = adj[v]
            if not N:
                continue
            # Fast path for k <= 1
            if k <= 1:
                cur_tag += 1
                tag_inN = cur_tag
                for u in N:
                    mark_inN[u] = tag_inN

                cur_tag += 1
                seen_tag = cur_tag
                comp_cnt = 0
                for u in N:
                    if mark_inN[u] != tag_inN or seen[u] == seen_tag:
                        continue
                    comp_cnt += 1
                    stack = [u]
                    seen[u] = seen_tag
                    while len(stack) > 0:
                        x = stack.pop()
                        neibours = adj[x]
                        for w in neibours:
                            if mark_inN[w] == tag_inN and seen[w] != seen_tag:
                                seen[w] = seen_tag
                                stack.append(w)
                tau[v] = comp_cnt
                continue

            # Mark neighborhood & compute local degrees
            cur_tag += 1
            tag_inN = cur_tag
            for u in N:
                mark_inN[u] = tag_inN
                removed[u] = 0
            for u in N:
                cnt = 0
                for w in adj[u]:
                    if mark_inN[w] == tag_inN:
                        cnt += 1
                deg_local[u] = cnt

            # Local k-core peeling
            peel_stack = [u for u in N if deg_local[u] < k]
            for u in peel_stack:
                removed[u] = 1
            while peel_stack:
                x = peel_stack.pop()
                for w in adj[x]:
                    if mark_inN[w] != tag_inN or removed[w]:
                        continue
                    deg_local[w] -= 1
                    if deg_local[w] == k - 1:
                        removed[w] = 1
                        peel_stack.append(w)

            # Count surviving components (all have local degree >= k)
            cur_tag += 1
            seen_tag = cur_tag
            comp_cnt = 0
            for u in N:
                if mark_inN[u] != tag_inN or removed[u] or seen[u] == seen_tag:
                    continue
                comp_cnt += 1
                stack = [u]
                seen[u] = seen_tag
                while len(stack) > 0:
                    x = stack.pop()
                    neighbours = adj[x]
                    for w in neighbours:
                        if (mark_inN[w] == tag_inN and not removed[w]
                                and seen[w] != seen_tag):
                            seen[w] = seen_tag
                            stack.append(w)

            tau[v] = comp_cnt

        return tau

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
