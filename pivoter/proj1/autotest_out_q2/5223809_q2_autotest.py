#!/usr/bin/env python3
# Auto-generated for 5223809

STUDENT_ID = "5223809"
STUDENT_NAME = "Yunchuan Li"

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

        # class My_Flat_PQ(object):
        #     def __init__(self):
        #         (self.d, self.b, self.D, self.p) = self.initialize_flat_array(G)
        #         self.pointer = 0

        def initialize_flat_array(G):
            n = G.vertex_num
            # initialize (d, b, D, p, G)
            
            # initialize array d
            # array d store deg of each vertices, the index of d is vertex id
            d = [len(G.adj_list[vid]) for vid in range(n)]
            max_deg = max(d)

            # bin sort d array by initialize n bins
            # array bin_D is a 2D array that store list of vertex id
            # each array in bin_D represent vertex with the same number of degree
            # the index of bin_D is the deg of vertex
            # add 1 for zero degree
            bin_D = [list() for _ in range(max_deg + 1)]
            for vid in range(len(d)):
                deg = d[vid]
                # insert vid into corresponding deg bin
                bin_D[deg].append(vid)

            # print(f"bin_D = {bin_D}")
            
            # Compute starting positions for each bin
            # initialize b array
            # array b count the 
            # Count how many vertices have each degree
            b = [0]
            total_len = 0
            x = 0
            for deg_bin in bin_D[:]:
                # print(f"there are {len(deg_bin)} nodes with deg{x}")
                x +=1
                total_len += len(deg_bin)
                b.append(total_len)

            # construct D array and p array according to bin_D
            D = []
            p = [0 for _ in range(n)]
            pos = 0
            for deg_bin in bin_D[:]:
                for vid in deg_bin:
                    D.append(vid)
                    # vid is at pos index in D
                    p[vid] = pos
                    pos += 1
            
            # print(f"d = {d}")
            # print(f"b = {b}")
            # print(f"D = {D}")
            # print(f"p = {p}")
            # exit()
            return (d, b, D, p)
            
        def core_decomposition(G):
            n = G.vertex_num
            (d, b, D, p) = initialize_flat_array(G)

            cn = [1 for _ in range(n)]
            # for each u in order do
            for i in range(n):
                v = D[i]
                cn[v] = d[v] # record
                # for each v in nbr(u) with d(v) > d(u) do
                for u in G.adj_list[v]:
                    if d[u] > d[v]:
                        du = d[u]   # degree of u
                        pu = p[u]   # position of u in D
                        pw = b[du]  # start position of degree of u in D
                        w = D[pw]   # vertex w which is at the start position of degree of u in D
                        if u != w:
                            # swap position of u and w in D
                            D[pu] = w
                            D[pw] = u
                            # update position of u and w in p for D
                            p[u] = pw
                            p[w] = pu
                        b[du] += 1
                        d[u] -= 1
            # print(cn)
            # print(f"dv\t - {[d[vid] for vid in D]}")
            # # print(f"D\t - {[G.id2label[v_id] for v_id in D]}")
            # print(f"D\t - {D}")
            # exit()
            return (d, b, D, p)
        
        class My_Sub_Graph(object):
            def __init__(self, u):
                self.u = u
                self.edge = []
                self.node = []
                self.node_label = []

                self.adj_list = []
                self.k_core = [-1 for _ in range(G.vertex_num)]

                self.label2id = {}
                self.id2label = []
                self.vertex_num = 0

            # O(1)
            def insertNode(self, v_label):
                self.node_label.append(v_label)

            # O(number of vertex)
            def constructNode(self):
                self.vertex_num = len(self.node_label)
                if self.vertex_num > 0:
                    self.node = [i for i in range(self.vertex_num)]
                    self.label2id = dict(zip(self.node_label, self.node))
                    self.id2label = dict(zip(self.node, self.node_label))

            # O(1)
            def insertEdge(self, v_label, w_label):
                # subGraph[u]
                # print(f"insert {self.u} {v_label}->{w_label}")
                vid = self.label2id[v_label]
                wid = self.label2id[w_label]
                self.edge.append([vid, wid])
            
            # O(number of edges)
            def constructGraph(self):
                self.adj_list = [list() for _ in range(self.vertex_num)]
                for u, v in self.edge:
                    self.adj_list[u].append(v)
                    self.adj_list[v].append(u)
        

        def TriangleCountInsert(G):
            n = G.vertex_num
            subGraph = [My_Sub_Graph(i) for i in range(n)]
            
            # 1. Orientation step, Construct Directed Graph for Triangle Count (directed graph stored as lists)
            # Time Complexity O(m)
            directed_adj = [[] for _ in range(n)]
            seenbitmap = [[False] * n for _ in range(n)]
            deg = [len(G.adj_list[i]) for i in range(n)]
            # for u_lab, u_id in G.vertex_dict.items():
            for u_id in range(n):
                # for v_lab in G.adj_list[u_id]:
                    # v_id = G.vertex_dict[v_lab]
                for v_id in G.adj_list[u_id]:
                    if deg[v_id] < deg[u_id] or (deg[v_id] == deg[u_id] and v_id > u_id):
                        # if v_id not in directed_adj[u_id]:          # avoid duplicates
                        #     directed_adj[u_id].append(v_id)
                        if seenbitmap[u_id][v_id] == False:
                            directed_adj[u_id].append(v_id)
                            seenbitmap[u_id][v_id] = True

            # 2.1 Triangle Counting Fetch Neighbour Nodes for all Nodes
            # O(∑(u, v) * 1) = O(m)
            for u in range(n):
                for v in directed_adj[u]:
                    subGraph[u].insertNode(v)
                    subGraph[v].insertNode(u)
            # what we construct here is the same as what we have just inserted
            # O(number of inserts) = O(m)
            for u in range(n):
                subGraph[u].constructNode()
            
            # 2.2 Triangle enumeration with a reusable bitmap
            # O(∑(u, v) * min(deg(u), deg(v))) = O(m^1.5)
            bitmap = [False] * n         # Boolean array
            modified = []                # list of ids whose bitmap entry became True
            triangles = 0
            for u in range(n):
                # mark neighbors of u
                for v in directed_adj[u]:
                    bitmap[v] = True
                    modified.append(v)
                # check directed wedges (u,v,w)
                for v in directed_adj[u]:
                    for w in directed_adj[v]:
                        if bitmap[w]:
                            triangles += 1
                            # print(G.id_to_raw[u], G.id_to_raw[v], G.id_to_raw[w])
                            # print(u, v, w)
                            # Triangle u=>v=>w=>u=>v, edge in subgraph is undirected after constructed
                            subGraph[u].insertEdge(v, w)
                            subGraph[v].insertEdge(w, u)
                            subGraph[w].insertEdge(u, v)

                # reset bitmap entries touched in this round
                for v in modified:
                    bitmap[v] = False
                modified.clear()
            return subGraph
        
        # BFS search to group the node ≥ k in the subGraph G
        # O(number of edge in subGraph)
        def BFS_joint(k, G, d, b, D, p):
            n = G.vertex_num
            visited = [False] * n
            cluster_list = []
            for uid in G.node:
                if d[uid] < k: continue
                if visited[uid]: continue
                cluster = []
                q = deque()
                q.append(uid)
                cluster.append(uid)
                visited[uid] = True
                
                while len(q) > 0:
                    sid = q.popleft()
                    for nbr_id in G.adj_list[sid]:
                        if d[nbr_id] < k: continue
                        if visited[nbr_id]: continue
                        q.append(nbr_id)
                        cluster.append(nbr_id)
                        visited[nbr_id] = True
                cluster_list.append(cluster)
            return cluster_list
        
        subGraph = []
        subGraph = TriangleCountInsert(G)
        # print(f"k = {1}")
        sd = []
        for uid in range(len(subGraph)):
            # if uid != 5: continue
            subG = subGraph[uid]
            # handle edge case where no edges or isolated vertices
            if len(subG.node) == 0 or len(subG.edge) == 0:
                sd.append(0)
                continue
            subG.constructGraph()
            # print(subG.node)
            # print([subG.id2label[vid] for vid in subG.node])
            # print(subG.adj_list)
            # print([[subG.id2label[nbr] for nbr in adj] for adj in subG.adj_list])
            (d, b, D, p) = core_decomposition(subG)
            # # calculate k-core for each subG.node
            # print(f"k_core d \t - {d}")
            # print(f"index d  \t - {[subG.id2label[i] for i in range(len(d))]}")
            # print(f"k_core D\t - {[d[vid] for vid in D]}")
            # print(f"index D \t - {[subG.id2label[v_id] for v_id in D]}")
            # exit()
            outputlist = BFS_joint(k, subG, d, b, D, p)
            
            # print(f"{uid}\t - {[[subG.id2label[nbr] for nbr in adj] for adj in outputlist]}")
            # print(f"\t - {[len([subG.id2label[nbr] for nbr in adj]) for adj in outputlist]}")
            sd.append(len(outputlist))
            # exit()
        # print(f"[kCore] {k} - {sd}")
        # exit()

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
