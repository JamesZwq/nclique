#!/usr/bin/env python3
# Auto-generated for 5447882

STUDENT_ID = "5447882"
STUDENT_NAME = "Yijia Hou"

# ======= 学生代码 =======
################################################################################
# You can import any Python Standard Library modules~
from collections import deque
################################################################################

class kCoreBaseStructuralDiversity(object):
    def __init__(self):
        pass

    @staticmethod
    def find_component(adj_list):
        # Use a list to store the connected components found during the query
        component_result = []
        # `visited` indicates which nodes have already been visited; these visited nodes will not be considered again as separate connected components
        visited = {}
        for node in adj_list:
            visited[node] = False
        # Traverse all nodes and skip those that have already been visited
        for node in adj_list:
            if visited[node] == True:
                continue
            # A node that hasn't been visited represents a new connected component, so create a list to store this component
            component = []
            visited[node] = True
            stack = [node]
            # Use a stack for iteration, continuously retrieve adjacent nodes and mark them as visited
            while stack:
                cur_node = stack.pop()
                component.append(cur_node)
                for adj_node in adj_list[cur_node]:
                    if visited[adj_node] == False:
                        stack.append(adj_node)
                        visited[adj_node] = True
            # Store each `component` into the result list
            component_result.append(sorted(component))
        return component_result

    @staticmethod
    def neighbor_induced_graph(adj_list,neighbors):
        # Create `result` to store the neighbor-induced subgraphs
        result={}
        neighbors_set=set(neighbors)
        # A neighbor-induced subgraph is the subgraph formed by the mutual connections among the neighbor nodes
        for neighbor in neighbors:
            # Traverse all neighbor nodes, and for each neighbor node A, find its neighbor B; if B is also a member of the neighbor set, add the edge between A and B to the neighbor-induced subgraph
            sub_neighbors=[]
            for adj_node in adj_list[neighbor]:
                if adj_node in neighbors_set:
                    sub_neighbors.append(adj_node)
            # Add the result to the `result` list
            result[neighbor]=sub_neighbors
        return result

    @staticmethod
    def update_subgraph(ni_subgraph,node,degrees):
        # Update the neighbor-induced subgraph of a node by deleting the node and all its adjacent edges, and update the degree of related nodes accordingly
        neighbors=ni_subgraph[node]
        if neighbors:
            for neighbor in ni_subgraph[node]:
                # Remove this node from the neighbor lists of all its neighboring nodes
                if neighbor in ni_subgraph:
                    ni_subgraph[neighbor].remove(node)
                    # Update the `degrees` dictionary
                    degrees[neighbor] -= 1
        # Directly remove this node from the neighbor-induced subgraph
        del ni_subgraph[node]
        del degrees[node]

        return ni_subgraph,degrees


    @staticmethod
    def compute_k_core(ni_subgraph,k):
        # The k-core is obtained by repeatedly removing nodes with degree < k and their adjacent edges, resulting in the final remaining subgraph
        degrees={}
        visited={}
        # Create the `visited` dictionary to prevent nodes from being visited multiple times
        for node in ni_subgraph:
            visited[node]=False
        # Initialize the `degrees` dictionary
        for node, neighbors in ni_subgraph.items():
            degrees[node]=len(neighbors)

        # Use a queue to find the initial batch of nodes with degree less than k and add them to the queue
        queue=deque()
        for node,degree in degrees.items():
            if degree<k:
                queue.append(node)
                visited[node]=True

        while queue:
            # As long as the queue contains nodes, pop them out and process the corresponding edges while updating the `degrees` dictionary
            cur_node = queue.popleft()
            neighbors = list(ni_subgraph[cur_node])
            # As nodes and edges are removed, the neighbor-induced subgraph gradually shrinks
            ni_subgraph,degrees=kCoreBaseStructuralDiversity.update_subgraph(ni_subgraph,cur_node,degrees)
            # Traverse again; as long as there are nodes with degree < k that haven't been visited, continue adding them to the queue to perform `update_subgraph`
            for neighbor_node in neighbors:
                if neighbor_node in degrees:
                    if degrees[neighbor_node]<k and visited[neighbor_node]==False:
                        queue.append(neighbor_node)
                        visited[neighbor_node]=True

        return ni_subgraph


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

        adj_list=G.adj_list

        # Compute τ_k(v) for each node v
        for vertex in range(n):
            neighbors=adj_list[vertex]
            # If there are no neighbor nodes, it means that τ_k(v) for this node must be 0
            if not neighbors:
                sd[vertex]=0
                continue

            # First compute the neighbor-induced subgraph
            ni_subgraph=kCoreBaseStructuralDiversity.neighbor_induced_graph(adj_list,neighbors)

            # Then compute the k-core of the neighbor-induced subgraph
            k_core = kCoreBaseStructuralDiversity.compute_k_core(ni_subgraph,k)

            # Perform connected component check
            result_k_core = kCoreBaseStructuralDiversity.find_component(k_core)

            # Add the result to the `sd` list
            sd[vertex] = len(result_k_core)

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
