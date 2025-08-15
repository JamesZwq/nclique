#!/usr/bin/env python3
# Auto-generated for 5507776

STUDENT_ID = "5507776"
STUDENT_NAME = "Shuyuan Xing"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque, defaultdict
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
        
    
        # find k + 1 core graph of graph G
        k_core_filter = kCoreBaseStructuralDiversity.find_k_core_graph(G.adj_list, k+1)
        
        # find all the components of the K + 1 core graph nodes
        total_comp, total_comps = kCoreBaseStructuralDiversity.find_num_components(G.adj_list, k_core_filter)
        
        # traverse all the comp for the k + 1 graph nodes
        for comp in total_comps:
            # traverse all the nodes in one comp
            for cur_node in comp:
                
                # find cur_node's neighbours
                nbr = set(G.adj_list[cur_node])
                
                # edge case: if cur_node have 0 neighbours, return 0
                if len(nbr) == 0:
                    sd[cur_node] = 0
                    continue
                
                # find cur_node's neighbor induced subgraph
                nbr_graph = kCoreBaseStructuralDiversity.find_nbrs(G.adj_list, cur_node)
                
                # find k core graph node in the cur_node's neighbor induced subgraph
                k_core_node = kCoreBaseStructuralDiversity.find_k_core_graph(nbr_graph, k)

                # compute the number of comps of k core graph node in the cur_node's neighbor induced subgraph
                num_comp, comps = kCoreBaseStructuralDiversity.find_num_components(nbr_graph, k_core_node)
                
                # store the number of comps
                sd[cur_node] = num_comp
        return sd


    ################################################################################
    # You can only define any auxiliary functions in class kCoreBaseStructuralDiversity
    ################################################################################
    @staticmethod
    def find_k_core_graph(adj_list, k):
    ###
    #   Function: K-core decomposition
    #   adj_list: Graph structure
    #   k: k-core number
    ###
    
        # Store deg of nodes, remove node and flag
        all_deg = defaultdict(int)
        remove_queue = deque()
        removed_flag = set() 
        
        # traverse graph whether graph structure is list or dict
        nodes = kCoreBaseStructuralDiversity.iterate_nodes(adj_list)
        
        # build all nodes degree and remove the node which degree is less than k
        for v in nodes:
            all_deg[v] = len(adj_list[v])
            if all_deg[v] < k:
                remove_queue.append(v)
                removed_flag.add(v)
        
        while remove_queue:
            # get cur remove node
            cur_node = remove_queue.popleft()
            
            # traverse cur node neighbours
            for u in adj_list[cur_node]:
                # if neighbour is not be removed
                if u not in removed_flag:
                    # degree - 1
                    all_deg[u] -= 1
                    # if neighbour degree less than k, delete neighbour
                    if all_deg[u] < k:
                        remove_queue.append(u)
                        removed_flag.add(u)
        
        # get K-core nodes   
        k_core_nodes = set(nodes) - removed_flag

        return k_core_nodes
    
    @staticmethod
    def find_nbrs(adj_list, node_id):
        ###
        #   Function: Get the neighbor induced subgraph
        #   adj_list: Graph structure
        #   node_id: target node
        ###
        
        # Initial nbr induced graph
        nbr_subgraph = defaultdict(list)
        
        # get node_id's neighbour
        nbrs = set(adj_list[node_id])
        
        # traverse node_id's neighbour
        for v in nbrs:
            v_nbrs = []
            # find neighbours which have edges between them
            for u in adj_list[v]:
                if u in nbrs:
                    v_nbrs.append(u)
            nbr_subgraph[v].extend(v_nbrs)
            
        return nbr_subgraph
    
    @staticmethod
    def find_num_components(adj_dict, node_set):
        ###
        #   Function: compute components
        #   adj_dict: Graph structure
        #   node_set: node list to find number of components
        ###
        
        # edge case for none node list
        if not node_set:
            return 0, []
        
        # Initial visited list, components code, number of components
        visited = set()
        num_components = 0
        components_node = []
        
        # traverse node list
        for v in node_set:
            # only visited node which is not be visited
            if v not in visited:
                visited.add(v)
                
                # found a new component
                num_components += 1
                # add first node to queue and component node list
                cur_queue = deque([v])
                cur_components = {v}
                
                # find others node in the current component
                while cur_queue:
                    # get next node
                    u = cur_queue.popleft()
                    
                    # traverse node u's neighbours
                    for i in adj_dict[u]:
                        
                        # find next node in the current component
                        if i not in visited and i in node_set:
                            # add it to visited list, queue and component node list
                            visited.add(i)
                            cur_queue.append(i)
                            cur_components.add(i)
                # add all the nodes in the current component
                components_node.append(cur_components)
        return num_components, components_node
    
    
    @staticmethod
    def iterate_nodes(adj):
        ###
        #   Function: Graph structure type identifier
        #   adj: Graph structure
        ###
        if isinstance(adj, dict):
            return adj.keys()          
        elif isinstance(adj, list):
            return range(len(adj))

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
